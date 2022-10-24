"""
    create_run(run_name::AbstractString)

Generates a generic input file ("run_name.in") loaded with options for the run in the current
directory.
"""
function create_run(run_name::AbstractString)
    # I would strongly recommend requiring the user to specify the
    # file extension rather than assuming it's *.in. -BF
    open(string(run_name, ".in"), write=true) do io
        lns = [
            "RUN_TYPE 1",
            "SITES 1:36",
            "SITE_LIST sitelist.txt",
            "AO 1",
            "RADIUS 2.2",
            "CONTINUE_FROM lastrun.mat"
        ]
        for s in lns
            println(io, s)
        end
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
    import_VASP(d::AbstractString="") -> fermi::Float64, geo::Crystal{3}, kpt::KPointGrid{3}, super::AtomList{3}, wave::ReciprocalWaveFunction{3,Float32}

Searches in the specified directory (default is the current directory) for VASP files
OUTCAR, POSCAR, KPOINTS, and WAVECAR, and extracts relevant information.
"""
function import_VASP(d::AbstractString="")
    fermi = get_fermi(string(d,"OUTCAR"))
    geo = readPOSCAR(string(d,"POSCAR"))
    kpt = readKPOINTS(string(d,"KPOINTS"))
    wave = readWAVECAR(string(d,"WAVECAR"))
    # Creates an AtomList{3} supercell from POSCAR and KPOINTS
    super = supercell(geo.atoms,kpt.grid)
    return (fermi, geo, kpt, super, wave)
end

"""
    read_GCOEFF(emin::Real, emax::Real) -> (g_indices::Vector{Int64}, occupied_coefficients::Matrix{ComplexF64}, kptlist::Vector{Vector{Float64}})

Reads a GCOEFF.txt file within the specified energy range. Returns the unique G vectors, complex coefficients, and k-point list corresponding to occupied states.
"""
function read_GCOEFF(emin::Real, emax::Real)
    open("GCOEFF.txt","r") do io
        num_spin_states = parse(Int,readline(io))
        num_kpt = parse(Int,readline(io))
        num_band = parse(Int,readline(io))
        # skips the next 6 lines, which correspond to a, b, c, a*, b*, c*
        for i in 1:6
            readline(io);
        end
        # Regarding spin states
        if num_spin_states == 2
            num_spin_up = 0
            num_spin_down = 0
        end
        # Initialize containers
        num_occ_states = 0
        kptlist = Array{Vector{Float64}}(undef,0)
        #current_G = Array{Vector{Int}}(undef,0)
        occ_coeff = []
        # Loop through spin states, kpts, band
        for s in 1:num_spin_states
            println("Reading spin state ", s)
            for k in 1:num_kpt
                current_kpt = parse.(Float64,split(readline(io)))
                # Counter for occupied states in kpt
                occ_states_k = 0
                for b in 1:num_band
                    current_band,num_pw = parse.(Int,split(readline(io)))
                    # Read in Energy, Filling
                    ln = split(readline(io))
                    # Need if statement somewhere to check if we've reached the end of the bands
                    energy, filling = parse.(Float64, [ln[2],ln[6]])
                    # Reads in info for only occupied states
                    if energy >= emin && energy <= emax
                        # Counts spin up/down states
                        if num_spin_states == 2
                            if s == 1
                                num_spin_up += 1
                            else
                                num_spin_down += 1
                            end
                        end
                        # Counts occupied states, updates kptlist
                        occ_states_k = occ_states_k + 1
                        num_occ_states = num_occ_states + 1
                        push!(kptlist,current_kpt)
                        # Loops through pw in each band and stores it in a array of tuples
                        # where the NamedTuple is (state = Int, G = Vector{Int}, coeff = ComplexF64)
                        # This probably needs a better method here.
                        for p in 1:num_pw
                            ln = split(readline(io))
                            temp_g = parse.(Int,ln[1:3])
                            temp_occ_coeff = complex(parse(Float64,ln[5]),parse(Float64,ln[7]))
                            push!(occ_coeff,(state = num_occ_states, G = temp_g, coeff = temp_occ_coeff))
                        end
                    else
                        for p in 1:num_pw
                            readline(io)
                        end
                    end
                end
                println("For K-point # ", k, ", Spin State ", s, ", found ", occ_states_k, " Occupied States.")
            end
        end
        println("Total occupied states: ", num_occ_states,)
        # Post processing of occ_coeff. Finds the unique G
        g_indices = unique(i->occ_coeff[i].G,eachindex(occ_coeff))
        unique_G = fill!(Vector{Vector{Int64}}(undef,length(g_indices)),[0,0,0])
        occupied_coefficients = fill!(Array{ComplexF64}(undef,length(g_indices),num_occ_states),0)
        # Finds unique G in occ_coeff, then fills an array at index [G, number_of_the_occupied_state] with the corresponding coefficient.
        # Quite inefficient. This will need optimization, but for now it works.
        for i in eachindex(g_indices)
            indices = findall(x->x==occ_coeff[g_indices[i]].G, [occ_coeff[j].G for j in eachindex(occ_coeff)])
            for j in indices
                occupied_coefficients[i,occ_coeff[j].state] = occ_coeff[j].coeff
                unique_G[i] = occ_coeff[j].G
            end
        end
        return (unique_G, occupied_coefficients, kptlist)
    end
end

# Repeated kpoints marked with 0; unrepeated marked with 1.
function track_kpoint_repeating(kptlist)#::KPointList{3})
    kpoint_repeating = zeros(Bool,length(kptlist))
    kpoint_repeating[1] = 1
    for i in 2:length(kptlist)
        if kptlist[i] != kptlist[i-1]
            kpoint_repeating[i] = 1
        end
    end
    return kpoint_repeating
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
        orb_param = OrbitalParams(parse(Int,current[2]), parse(Int,current[3]), get(orbs,current[6],3), parse(Int,current[5]), parse(Float64,current[7]), parse(Float64,current[8]), parse(Float64,current[9]), parse(Float64,current[10]), parse(Float64,current[11]))
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
        # Removes first column of atomic labels, if any
        has_label = 0
        if length(ln[1])==4
            has_label = 1
        end
        for i in ln
            push!(sitelist,parse.(Float64,i[1+has_label:3+has_label]))
        end
    end
    return sitelist
end