"""
    import_VASP(d::AbstractString="")
        -> Tuple{
            NamedTuple{(:fermi, :alphabeta),NTuple{2,Float64}},
            Crystal{3},
            PlanewaveWavefunction{3,ComplexF32}
        }

Searches in the specified directory (defaults to the current directory if none is given) for VASP
`OUTCAR`, `POSCAR`, `KPOINTS`, and `WAVECAR` files, returning a `Tuple` with all of the data needed
for a DFT-raMO run.
"""
function import_VASP(directory::AbstractString="")
    (isdir(directory) || isempty(directory)) || error("$directory is not a directory.")
    fermi = get_fermi(directory)
    geo = readPOSCAR(directory)
    wave = readWAVECAR(directory, quiet = true)
    kpt = parse.(Int, split(readlines(string(directory, "KPOINTS"))[4]))
    # Use a Crystal to lazily reference the supercell
    xtal = set_transform!(Crystal(geo), kpt)
    xtal = PeriodicAtomList(xtal)
    return (fermi=fermi, xtal=xtal, geo=geo, kpt=kpt, wave=wave)
end

"""
    import_abinit(file)

Reads a abinit binary wavefunction output (usually ending in `_WFK`) and gets the Fermi energy, the
crystal geometry, the k-point mesh, and the wavefunction data.
"""
function import_abinit(io::IO)
    h = Electrum.read_abinit_header(io)
    seekstart(io)
    isdiag(h[:kpt]) || error("Currently, non-diagonal k-point lattices are unsupported.")
    return (
        fermi = h[:fermi],
        xtal = Crystal(h),
        geo = PeriodicAtomList(Crystal(h)),
        kpt = diag(h[:kptrlatt]),
        wave = read_abinit_WFK(io)["wavefunction"]
    )
end

import_abinit(file) = open(import_abinit, file)

"""
    get_occupied_states(wave::PlanewaveWavefunction, energy::Real) ->
        occ_state_array::Array{DFTraMO.OccupiedState}

Returns an array of OccupiedState where the energy of the wavefunction is less
than the specified energy (usually the fermi energy). The OccupiedState array
has the dimensions n (number of occupied states) by m (number of nonzero
planewaves). Each OccupiedState element holds information about the coefficient,
kpoint, and G vector.
"""
function get_occupied_states(wave::PlanewaveWavefunction, emin::Real, emax::Real)
    hkl_list = [SVector(v.I) for v in FFTBins(wave)]

    # Filters occupied states within energy range
    num_occ_states = emin .< wave.energies .< emax # Bit array
    occ_states = reshape(
        [
            (wave.data[n], wave.kpoints.points[CartesianIndices(wave.data)[n].I[3]].point, hkl_list[CartesianIndices(wave.data)[n].I[1]])
            for n in eachindex(wave.data)
        ],
        size(wave.data)
    )
    occ_states = [occ_states[n,:,:,:][num_occ_states] for n in 1:size(occ_states)[1]]
    
    # Filters planewaves that are nonzero throughout all kpoints and bands
    nonzero_coeff = [iszero(occ_states[n][m][1]) for n in eachindex(occ_states), m in eachindex(occ_states[1])]
    occ_pw = [sum(nonzero_coeff[n,:]) != sum(num_occ_states) for n in 1:size(nonzero_coeff)[1]]

    # Create an n occupied states by m occupied planewaves array that stores the
    # coefficient, kpoint, and corresponding G vector.
    occ_states = occ_states[occ_pw]
    occ_states = mapreduce(permutedims, vcat, [occ_states[n] for n in eachindex(occ_states)])

    # coeff is a matrix of tuples with dimensions occ_planewave x occ_states
    # (occ_coeff, kpt, hkl)
    coeff = Matrix{ComplexF32}(undef, size(occ_states))
    kpt = Array{SVector{3, Float64}}(undef, size(occ_states))
    hkl_list = Array{SVector{3, Int64}}(undef, size(occ_states))
    for n in eachindex(occ_states)
        coeff[n] = occ_states[n][1]
        kptvect = Vector{Float64}(undef,3)
        # correct for close-to-zero kpt vectors
        for k in eachindex(occ_states[n][2])
            if abs(occ_states[n][2][k]) < 0.000001
                kptvect[k] = 0
            else
                kptvect[k] = occ_states[n][2][k]
            end
        end
        kpt[n] = SVector{3, Float64}(kptvect)
        hkl_list[n] = occ_states[n][3]
    end

    return OccupiedStates(coeff, kpt, hkl_list)
end

get_occupied_states(x::raMODFTData, emin::Real, emax::Real) = get_occupied_states(x.wave, emin, emax)
get_occupied_states(r::raMORuns) = get_occupied_states(PlanewaveWavefunction(r), r.emin, r.emax)

"""
    read_eht_params(paramsfile::AbstractString) -> mat::ehtParams 

Reads an eHtuner parameter file and returns it as an ehtParams object.
If argument is left empty, it will read the DFT_raMO_eht_parms.dat file by default.
"""
# Remove from export list
function read_eht_params(paramsfile::AbstractString="testfiles/DFT_raMO_eht_parms.dat")
    # Skips header
    ln = readlines(open(paramsfile,"r"))
    ln = ln[4:length(ln)]
    ln = filter(!isempty,split.(ln))
    unique_atoms_index = unique(i->mapreduce(permutedims,vcat,ln)[:,1][i],eachindex(mapreduce(permutedims,vcat,ln)[:,1]))
    mat = ehtParams(Matrix{OrbitalParams}(undef,length(ln),4))
    # Prefill with zeros
    for i in 1:length(mat.data)
        mat.data[i] = OrbitalParams(0,0,0,0,0.,0.,0.,0.,0.)
    end
    atom_counter = 0
    # Runs through list of unique atoms
    for i in eachindex(ln)
        current = ln[i]
        # Assigns atom_num,valence,l_quant,n_quant,IP,exp1,exp2,coeff1,coeff2
        orbs = Dict("s"=>0,"p"=>1,"d"=>2,"f"=>3)
        to_ang = 0.52917721092 # constant from converting from a.u. to angstrom
        orb_param = OrbitalParams(
            parse(Int,current[2]), # Atomic number
            parse(Int,current[3]), # Valence
            get(orbs,current[6],3), # l_quant
            parse(Int,current[5]), # n_quant
            parse(Float64,current[7]), # IP
            parse(Float64,current[8])/to_ang, # exp1
            parse(Float64,current[9])/to_ang, # exp2
            parse(Float64,current[10]), # coeff1
            parse(Float64,current[11]) # coeff2
            )
        if i in unique_atoms_index
            atom_counter = atom_counter + 1
        end
        # Info stored in matrix of (unique atoms, orbital)
        mat.data[atom_counter, get(orbs,current[6],3)+1] = orb_param
    end
    return mat
end

"""
    read_site_list(filename::AbstractString) -> site_list::Vector{Vector{Float64}}

Reads in a txt file with coordinates, typically for specifing midpoints for isolobal bonds or cage states.
"""
#TODO: predict length of list, change to matrices instead of vector of vector (will need to check other method)
function read_site_list(filename::AbstractString)
    sitelist = Vector{Vector{Float64}}(undef,0)
    open(filename,"r") do io
        ln = readlines(io)
        ln = filter(!isempty,split.(ln))
        length(ln[1]) > 4 && error("Check that your file is in the format 'Atom  0.0  0.0  0.0'")
        # Removes first column of atomic labels, if any
        has_label = (length(ln[1]) == 4)
        for i in ln
            push!(sitelist,parse.(Float64,i[1+has_label:3+has_label]))
        end
    end
    return sitelist
end

"""
    import_checkpoint(filename::AbstractString)
        -> (
            psi_previous::Array{ComplexF32},
            num_electrons_left::Int,
            num_raMO::Int
            )

Imports a matrix of remainder coefficients to use for raMO runs, as well as the number of remaining
electrons and the number of the raMO in the sequence.
"""
function import_checkpoint(filename::AbstractString)
    num_electrons_left = parse(Int,split(filename, ['.', '_'])[end-1])
    num_raMO = parse(Int,split(filename, ['.', '_'])[end-2])
    psi_previous = open(filename, "r") do io
        dims = [read(io, Int64) for n in 1:3] # first three Ints are the dimensions of the matrix
        psi_previous = Array{ComplexF32}(undef, (dims[1], dims[2], dims[3]))
        read!(io, psi_previous)
        return psi_previous
    end
    return (psi_previous, num_electrons_left, num_raMO)
end

"""
    import_raMO(filename::AbstractString) -> psi::Vector{ComplexF32}

Imports a vector of coefficients corresponding to each raMO.
"""
function import_raMO(filename::AbstractString)
    psi = open(filename, "r") do io
        sz = Int(stat(io).size/sizeof(ComplexF32))
        psi = Vector{ComplexF32}(undef, sz)
        read!(io, psi)
        return psi
    end
    return psi
end
