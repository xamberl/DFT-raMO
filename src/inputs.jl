"""
    read_run_yaml(run_name::AbstractString)

Loads options for the run specified by a file.
"""
function read_run_yaml(file::AbstractString, software::AbstractString="vasp")
    data = YAML.load_file(file)
    software == "vasp" ? dftinfo = import_VASP() : error("Other DFT packages are not implement yet.")

    # crayons for color
    cr_b = crayon"light_cyan"
    cr_y = crayon"yellow"
    cr_d = crayon"default"

    checkpoint = get(data, "checkpoint", nothing)
    if isnothing(checkpoint)
        println("No checkpoint file specified. Run will start from beginning conditions.")
    elseif !isfile(checkpoint)
        error(checkpoint, " does not exist. Check path to file.") : nothing
    else
        println("Using checkpoint file ", cr_b, checkpoint, cr_d)
    end 
    
    auto_psphere = get(data, "auto_psphere", false)
    typeof(auto_psphere) != Bool ? error("Invalid entry for auto_psphere.") : nothing
    println("Auto-Psphere: ", auto_psphere)
    
    runs = get(data, "runs", nothing)
    println("Number of runs: ", length(runs))
    runlist = Vector{RunInfo}(undef, length(runs))
    for n in eachindex(runs)
        println("Run ", n, ":")
        # Checks type, site_file, sites, radius, rsphere
        name = get(runs[n], "name", string("run_", n))
        println("   name: ", cr_b, name, cr_d)

        type = lowercase(get(runs[n], "type", ""))
        isempty(type) ? error("type for run ", n, " is empty") : nothing
        !in(type, union(keys(AO_RUNS), CAGE_RUNS, ["lcao"])) ? error("type ", type, " is an invalid entry.") : nothing
        println("   type: ", cr_b, type, cr_d)

        site_file = get(runs[n], "site_file", nothing)
        if isnothing(site_file)
            site_file = ""
        else
            # site_file is not needed for an AO run
            if in(type, keys(AO_RUNS))
                println("   AO type run. site_file will be ignored.")
                site_file = ""
            else
                !isfile(site_file) ? error(site_file, " does not exist. Check filename.") : nothing 
                if in(type, CAGE_RUNS)
                    site_list = read_site_list(site_file)
                elseif type == "salc"
                    salc_yaml = YAML.load_file(site_file)
                    site_list = get(salc_yaml, "salcs", nothing)
                end
                println("   site_file: ", cr_b, site_file, cr_d)
            end
        end

        sites = get(runs[n], "sites", nothing)
        isnothing(sites) ? error("sites cannot be empty.") : nothing
        site_final = Vector{Int}(undef, 0)
        if sites == "all"
            in(type, keys(AO_RUNS)) ? error("'all' is not valid for atomic orbital type runs.") : nothing
            sites_final = collect(1:length(site_list))
        else
            sites = split(sites, [' ', ','], keepempty=false)
            # Different requirements for different type of runs
            # If CAGE_RUNS, no element needs to be specified
            # If AO_RUNS, an element must be the first item in the list.
            if in(type, keys(AO_RUNS))
                e = string(sites[1])
                !in(e, Electrum.ELEMENTS) ? error(e, " is not a valid element.") : nothing
                valid_atoms = findall(x -> Electrum.name(x) == e, dftinfo.xtal.atoms)
                isempty(valid_atoms) ? error("Element ", e, "was not found in the system.") : nothing
                length(sites) > 1 ? sites_final = valid_atoms[parse_sites(sites[2:end])] : sites_final = valid_atoms
            else
                sites_final = parse_sites(sites)
                length(sites_final) > length(site_list) ? error("Number of sites specified exceeds site_list length.") : nothing
            end
        end
        println("   sites: ", cr_b, get(runs[n], "sites", nothing), cr_d)
        
        radius = get(runs[n], "radius", 0.0)
        if isnothing(radius)
            radius = 0.0
        else
            if in(type, keys(AO_RUNS))
                println("   Radius is ignored for atomic orbital type runs.")
            else
                typeof(radius) == String ? error("Radius must be a number.") : nothing
                println("   radius: ", cr_b, radius, " Å", cr_d)
            end
        end

        rsphere = get(runs[n], "rsphere", nothing)
        if isnothing(rsphere) 
            rsphere = 3.0
            println(cr_y, "   rsphere was not specified. Default value of 3.0 Å is applied.", cr_d)
        end
        typeof(rsphere) == String ? error("rsphere must be a number.") : Float64(rsphere)
        println("   rsphere: ", cr_b, rsphere, " Å", cr_d)

        runlist[n] = RunInfo(name, type, site_file, sites_final, radius, rsphere)
    end
    return(runlist, checkpoint, auto_psphere, dftinfo)
end

"""
    parse_sites(sites::Vector{SubString{String}}) -> site_final::Vector{Int}

Parses a portion of the "sites" lines in the yaml file to return a Vector{Int} with valid
indices for targets. e.g. ["1:3", "3", "18:2:20"] -> [1, 2, 3, 18, 20, 22]
"""
function parse_sites(sites::Vector{SubString{String}})
    site_final = Vector{Int}(undef, 0)
    for s in sites
        # process ranges of Ints, i.e. "1:3"
        if occursin(":", s)
            ln = parse.(Int, split(s, ':'))
            length(ln) == 2 ? site_final = vcat(site_final, collect(ln[1]:ln[2])) : nothing
            if length(ln) == 3
                check_range = collect(ln[1]:ln[2]:ln[3])
                typeof(check_range) == Vector{Float64} ? error("Atom range ", ln, "is invalid.") : nothing
                site_final = vcat(site_final, collect(ln[1]:ln[2]:ln[3]))
            end            
            # account for erroneous cases
            length(ln) == 1 || length(ln) > 3 ? error("Error in sites input: ", s, ". Please fix.") : nothing
        else
            site_final = vcat(site_final, parse(Int, s))
        end
    end
    return unique(site_final)
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
    psi = Vector{ComplexF64}(undef, size(file)[1])
    for i in eachindex(psi)
       real = file[i,1]
       file[i,2] == "+" ? n = 1 : n = -1
       imag = parse.(Float64,split(file[i,3],[' ', 'f', 'i', 'm'])[1])*n
       psi[i] = real+imag*im
    end
    return (psi)
end
