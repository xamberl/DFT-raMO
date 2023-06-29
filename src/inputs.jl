"""
    create_run(run_name::AbstractString)

Generates a generic input file ("run_name.in") loaded with options for the run in the current directory.
"""
function create_run(run_name::AbstractString)
    open(string(run_name,".in"),"w") do io
        println(io,
        "RUN_TYPE 1\nSITES 1:36\nSITE_LIST sitelist.txt\nAO 1\nRADIUS 2.2\nCONTINUE_FROM lastrun.mat")
    end
end

"""
    read_run_in(run_name::AbstractString)

Reads and loads options for the run specified by a file.
"""
function read_run_in(run_name::AbstractString)
    # initialize outputs
    run_type = 0
    sites = Vector{Int}(undef,0)
    sites_list = ""
    AO = 0
    radius = 0.0
    continue_from = ""
    # Begin parsing
    ln = readlines(run_name)
    for i in eachindex(ln)
        ln_temp = lstrip(ln[i])
        startswith(ln_temp, "RUN_TYPE") ? run_type = parse(Int,split( ln_temp)[2]) : nothing
        # Reads and parses SITES
        if startswith(ln_temp, "SITES")
            # Removes spaces and commas after SITE_LIST
            temp_sites = split(ln_temp,[' ',','])[2:length(split(ln_temp,[' ',',']))]
            temp_sites = filter(!isempty, temp_sites)
            # Checks for any site ranges
            for j in temp_sites
                temp = parse.(Int,split(j,':'))
                if length(temp) == 2
                    sites = vcat(sites,vcat(temp[1]:temp[2]))
                else
                    sites = vcat(sites,temp[1])
                end
            end
        end
        startswith(ln_temp, "SITE_LIST") ? sites_list = split( ln_temp)[2] : nothing
        startswith(ln_temp, "AO") ? AO = parse(Int,split( ln_temp)[2]) : nothing
        startswith(ln_temp, "RADIUS") ? radius = parse(Float64,split( ln_temp)[2]) : nothing
        startswith(ln_temp, "CONTINUE_FROM") ? continue_from = split( ln_temp)[2] : nothing
    end
    # Checks and warnings
    run_type == 0 ? @error("RUN_TYPE not found!") : nothing
    isempty(sites) ? @error("SITES not found!") : nothing
    isempty(sites_list) && run_type > 1 ? @error("SITE_LIST not found for hybrid run!") : nothing
    radius <= 0 && run_type > 1 ? @error("RADIUS not specified or is negative for hybrid run!") : nothing
    isempty(continue_from) ? @info("Starting afresh at the maximum valence electron count.") : @info("Continuing from ", continue_from)
    return run_type, sites, sites_list, AO, radius, continue_from
end

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
    wave = readWAVECAR(directory)
    kpt = parse.(Int, split(readlines(string(directory, "KPOINTS"))[4]))
    # Use a Crystal to lazily reference the supercell
    xtal = set_transform!(Crystal(geo), kpt)
    xtal = PeriodicAtomList(xtal)
    return (fermi, xtal, geo, kpt, wave)
end


"""
    get_occupied_states(wave::PlanewaveWavefunction, energy::Real) ->
        occ_state_array::Array{DFTraMO.OccupiedState}

Returns an array of OccupiedState where the energy of the wavefunction is less
than the specified energy (usually the fermi energy). The OccupiedState array
has the dimensions n (number of occupied states) by m (number of nonzero
planewaves). Each OccupiedState element holds information about the coefficient,
kpoint, and G vector.
"""
function get_occupied_states(wave::PlanewaveWavefunction, energy::Real)
    hkl_list = [SVector(v.I) for v in FFTBins(wave)]

    # Filters occupied states below specified energy
    num_occ_states = wave.energies .< energy # Bit array
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
    coeff = Array{ComplexF32}(undef, size(occ_states))
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
        kpt[n] = SVector{3, Float64}(kptvect) #occ_states[n][2]
        hkl_list[n] = occ_states[n][3]
    end

    return OccupiedStates(coeff, kpt, hkl_list)
end

"""
    read_eht_params(paramsfile::AbstractString) -> mat::ehtParams 

Reads an eHtuner parameter file and returns it as an ehtParams object.
If argument is left empty, it will read the DFT_raMO_eht_parms.dat file by default.
"""
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
function read_site_list(filename::AbstractString)
    sitelist = Vector{Vector{Float64}}(undef,0)
    open(filename,"r") do io
        ln = readlines(io)
        ln = filter(!isempty,split.(ln))
        length(ln[1]) > 4 ? error("Check that your file is in the format 'Atom  0.0  0.0  0.0'") : nothing
        # Removes first column of atomic labels, if any
        length(ln[1])==4 ? has_label = 1 : has_label = 0
        for i in ln
            push!(sitelist,parse.(Float64,i[1+has_label:3+has_label]))
        end
    end
    return sitelist
end

"""
    import_psi_previous(filename::AbstractString)

Imports a remainder matrix of coefficients to use for raMO runs.
"""
function import_psi_previous(filename::AbstractString)
    file = readdlm(filename, '\t')
    psi_previous = Array{ComplexF32}(undef, size(file))
    for i in eachindex(file)
       f = split(file[i],[' ', 'f', 'i', 'm'])
       real = parse(Float32,f[1])*10^parse(Float32,f[2])
       imag = parse(Float32,string(f[3],f[4]))*10^parse(Float32,f[5])
       psi_previous[i] = real+imag*im
    end
    num_electrons_left = parse(Int,split(filename, ['.', '_'])[end-1])
    num_raMO = parse(Int,split(filename, ['.', '_'])[end-2])
    return (psi_previous, num_electrons_left, num_raMO)
end

function import_psi(filename::AbstractString)
    file = readdlm(filename)
    psi = Vector{ComplexF32}(undef, size(file)[1])
    for i in eachindex(psi)
       real = file[i,1]
       file[i,2] == "+" ? n = 1 : n = -1
       imag = parse.(Float32,split(file[i,3],[' ', 'f', 'i', 'm'])[1])
       psi[i] = real+imag*im
    end
    return (psi)
end