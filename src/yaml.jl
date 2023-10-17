"""
    read_run_yaml(run_name::AbstractString)

Loads options for the run specified by a file.
"""
function read_run_yaml(file::AbstractString, software::AbstractString="vasp")
    data = YAML.load_file(file)
    software == "vasp" ? dftinfo = import_VASP() : error("Other DFT packages are not implement yet.")

    # crayons for color
    cr_b = crayon"light_cyan"
    cr_d = crayon"default"

    checkpoint = get(data, "checkpoint", nothing)
    if isnothing(checkpoint)
        println("No checkpoint file specified. Run will start from beginning conditions.")
    elseif !isfile(checkpoint)
        error(checkpoint, " does not exist. Check path to file.")
    else
        println("Using checkpoint file ", cr_b, checkpoint, cr_d)
    end 
    
    auto_psphere = get(data, "auto_psphere", false)
    typeof(auto_psphere) != Bool && error("Invalid entry for auto_psphere.")
    println("Auto-Psphere: ", auto_psphere)

    # Get energy range in Ha
    emin = get(data, "emin", nothing)
    if isnothing(emin)
        println("No minimum energy specified. All bands below emax will be included.")
        emin = minimum(dftinfo.wave.energies)
    else
        emin = parse_energy(software, emin)
    end
    emin > maximum(dftinfo.wave.energies) && error("emin cannot be greater than ", maximum(dftinfo.wave.energies), "Ha.")   
    emax = get(data, "emax", nothing)
    if isnothing(emax)
        println("No maximum energy specified. emax will be set to the Fermi energy.")
        emax = dftinfo.fermi.fermi
        software == "vasp" ? emax = emax*Electrum.EV2HARTREE : nothing
    else
        emax = parse_energy(software, emax)
    end
    emax < emin && error("emax cannot be less than emin.")
    println("Energy range: ", @sprintf("%.3f", emin), " to ", @sprintf("%.3f", emax), " Ha (", @sprintf("%.3f", emin*Electrum.HARTREE2EV), " to ", @sprintf("%.3f", emax*Electrum.HARTREE2EV), " eV)")
    
    runs = get(data, "runs", nothing)
    println("Number of runs: ", length(runs))
    runlist = parse_runs(runs, dftinfo)

    return(runlist, checkpoint, auto_psphere, dftinfo, emin, emax)
end

"""
    parse_runs(runs::Vector{Dict{Any, Any}}) -> runlist::Vector{RunInfo}

Parse the run section of the input yaml file and returns Vector{RunInfo}.
"""
function parse_runs(runs::Vector{Dict{Any, Any}}, dftinfo)
    cr_b = crayon"light_cyan"
    cr_y = crayon"yellow"
    cr_d = crayon"default"
    runlist = Vector{RunInfo}(undef, length(runs))
    for n in eachindex(runs)
        println("Run ", n, ":")
        # Checks type, site_file, sites, radius, rsphere
        name = get(runs[n], "name", string("run_", n))
        println("   name: ", cr_b, name, cr_d)

        type = lowercase(get(runs[n], "type", ""))
        isempty(type) && error("type for run ", n, " is empty")
        !in(type, union(keys(AO_RUNS), CAGE_RUNS, ["lcao"])) && error("type ", type, " is an invalid entry.")
        println("   type: ", cr_b, type, cr_d)

        site_file = get(runs[n], "site_file", "")

        # site_file is not needed for an AO run
        if in(type, keys(AO_RUNS))
            println("   AO type run. site_file will be ignored.")
            site_file = ""
        else
            !isfile(site_file) && error(site_file, " does not exist. Check filename.")
            if in(type, CAGE_RUNS)
                site_list = read_site_list(site_file) #
            elseif type == "lcao"
                lcao_yaml = YAML.load_file(site_file)
                site_list = get(lcao_yaml, "lcao", nothing)
            end
            println("   site_file: ", cr_b, site_file, cr_d)
        end

        sites = get(runs[n], "sites", nothing)
        isnothing(sites) ? error("sites cannot be empty.") : nothing
        if sites == "all"
            in(type, keys(AO_RUNS)) && error("'all' is not valid for atomic orbital type runs.")
            sites_final = collect(1:length(site_list))
        else
            sites = split(sites, [' ', ','], keepempty=false)
            # Different requirements for different type of runs
            # If CAGE_RUNS, no element needs to be specified
            # If AO_RUNS, an element must be the first item in the list.
            if in(type, keys(AO_RUNS))
                e = string(sites[1])
                !in(e, Electrum.ELEMENTS) && error(e, " is not a valid element.")
                valid_atoms = findall(x -> Electrum.name(x) == e, dftinfo.xtal.atoms)
                isempty(valid_atoms) && error("Element ", e, "was not found in the system.") 
                length(sites) > 1 ? sites_final = valid_atoms[parse_sites(sites[2:end])] : sites_final = valid_atoms
            else
                sites_final = parse_sites(sites)
                length(sites_final) > length(site_list) && error("Number of sites specified exceeds site_list length.")
            end
        end
        println("   sites: ", cr_b, get(runs[n], "sites", nothing), cr_d)
        
        radius = get(runs[n], "radius", 0.0)
        if isnothing(radius)
            radius = 0.0
        end
        if in(type, keys(AO_RUNS))
            println("   Radius is ignored for atomic orbital type runs.")
        else
            !isa(radius, Number) && error("Radius must be a number.")
            println("   radius: ", cr_b, radius, " Å", cr_d)
        end

        rsphere = get(runs[n], "rsphere") do
            println(cr_y, "   rsphere was not specified. Default value of 3.0 Å is applied.", cr_d)
            return 3.0
        end
        !isa(rsphere, Number) && error("rsphere must be a number.")
        println("   rsphere: ", cr_b, rsphere, " Å", cr_d)

        runlist[n] = RunInfo(name, type, site_file, sites_final, radius, rsphere)
    end
    return runlist
end

"""
    parse_energy(software::AbstractString, energy::AbstractString) -> n::Float64

Parse energy values of the input yaml files and returns energies in Ha.
"""
function parse_energy(software::AbstractString, energy::AbstractString)
    # Get energy range. Need to convert to Hartree
    ln = split(energy)
    length(ln) > 2 && error("Check energy range.")
    n = parse(Float64, ln[1])
    if lowercase(ln[2]) == "ev"
        n = n*Electrum.EV2HARTREE
    elseif !(lowercase(ln[2]) == "ha")
        error(ln[2], " is not a valid energy unit. Use Ha or eV.")
    end
    return n
end

"""
    parse_energy(software::AbstractString, energy::Real) -> n::Float64

Parse energy values of the input yaml files and returns energies in Ha.
"""
function parse_energy(software::AbstractString, energy::Real)
    if software == "vasp"
        e = Electrum.EV2HARTREE
        units = "eV"
    else
        e = 1
        units = "Ha"
    end
    warn("No energy units are specified. Defaulting to ", units, " (", software, ")")
    n = energy*e
    return n
end

"""
    parse_sites(sites::AbstractVector{<:AbstractString}) -> site_final::Vector{Int}

Parses a portion of the "sites" lines in the yaml file to return a Vector{Int} with valid
indices for targets. e.g. ["1:3", "3", "18:2:20"] -> [1, 2, 3, 18, 20, 22]
"""
function parse_sites(sites::AbstractVector{<:AbstractString})
    site_final = Vector{Int}(undef, 0)
    for s in sites
        # process ranges of Ints, i.e. "1:3"
        if occursin(":", s)
            ln = parse.(Int, split(s, ':'))
            length(ln) == 2 ? site_final = vcat(site_final, collect(ln[1]:ln[2])) : nothing
            if length(ln) == 3
                check_range = collect(ln[1]:ln[2]:ln[3])
                typeof(check_range) == Vector{Float64} && error("Atom range ", ln, "is invalid.")
                site_final = vcat(site_final, collect(ln[1]:ln[2]:ln[3]))
            end            
            # account for erroneous cases
            (length(ln) == 1 || length(ln) > 3) && error("Error in sites input: ", s, ". Please fix.")
        else
            site_final = vcat(site_final, parse(Int, s))
        end
    end
    return unique(site_final)
end