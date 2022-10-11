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
extract_VASP

Searches in the current directory for VASP files and extracts relevant information.
"""
function extract_VASP()
    fermi = get_fermi("OUTCAR")    
    geo = readPOSCAR("POSCAR")
    kpt = readKPOINTS("KPOINTS")
    wave = readWAVECAR("WAVECAR")
    # Creates an AtomList{3} supercell from POSCAR and KPOINTS
    super = supercell(geo.gen,kpt.grid)
    return fermi, geo, kpt, super, wave
end

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
            npw_temp = 0
            current_G = Array{Vector{Int}}(undef,0)
            current_occ_coeff = []
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
                            #temp_occ_coeff = Vector{Complex}(undef,0)
                            for p in 1:num_pw
                                # Counter to keep track of which G we're at
                                g_count = 1
                                ln = split(readline(io))
                                # If the G value is unique, append it to G, and temp_occ_coeff
                                if !(parse.(Int,ln[1:3]) in current_G)
                                    npw_temp = npw_temp+1
                                    push!(current_G, parse.(Int,ln[1:3]))
                                    push!(current_occ_coeff, complex(parse(Float64,ln[5]),parse(Float64,ln[7])))
                                # Otherwise, put it in the corresponding slot
                                else
                                    push!(current_occ_coeff, complex(parse(Float64,ln[5]),parse(Float64,ln[7])))
                                end
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
            println("Total occupied states: ", num_occ_states," npw_temp: ", npw_temp)
            return current_G, current_occ_coeff, kptlist
        end
    end

# Repeated kpoints marked with 0; unrepeated marked with 1.
function track_kpoint_repeating(kptlist::Vector{Vector{Float64}})
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
    read_eht_params(paramsfile::AbstractString) -> ehtParams 

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
    for i in 1:length(ln)
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