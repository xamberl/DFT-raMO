"""
create_run(run_name::AbstractString)

Generates a generic input file ("run_name.in") loaded with options for the run in the current directory.
"""
function create_run(run_name::AbstractString)
    open(string(run_name,".in"),"w") do io
        println(io,
        "run_type 1\nsites 1:36\nsite_list sitelist.txt\nAO 1\nradius 2.2\ncontinue_from lastrun.mat")
    end
end

"""
extract_VASP

Searches in the current directory for VASP files and extracts relevant information.
    """
    function extract_VASP()
        # Read Fermi energy from OUTCAR
        fermi = get_fermi("OUTCAR")
        
        # Get geometry from POSCAR and KPOINTS (WIP)
        geo = readPOSCAR("POSCAR")
        kpt = readKPOINTS("KPOINTS")
        # Creates an AtomList{3} supercell from POSCAR and KPOINTS
        super = supercell(geo.gen,kpt.grid)
        
        # Load in default number of atomic orbitals, number of electrons [WIP]
        
        # Prints info into a system.in file
        #==open("system.in","w") do io
            println(io,"fermi_energy ", fermi, "\nenergy_range -100:"," fermi")
        end==#
        return fermi, geo, kpt, super
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

# Creates overlap matrix with occupied coefficients
# If the kpoints are the same, they can overlap. Otherwise, they do not.
function make_overlap_mat(occ_coeff, kpoint_repeating::Vector{Bool})
    num_occ_states = length(kpoint_repeating)
    S = zeros(num_occ_states,num_occ_states)
    kpoint_zero = findall(x -> x == 0, kpoint_repeating)
    S[kpoint_zero,kpoint_zero] = occ_coeff[:,kpoint_zero]'*occ_coeff[:,kpoint_zero]
    return S
end