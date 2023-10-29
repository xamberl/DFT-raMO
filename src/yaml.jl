"""
    dftramo_run(filename::AbstractString, software::AbstractString="vasp")

Automatically runs DFTraMO from a configuration yaml file.
"""
function dftramo_run(filename::AbstractString, software::AbstractString="vasp")
    (runs, checkpoint, auto_psphere, dftinfo, emin, emax) = read_run_yaml(filename, software)
    occ_states = get_occupied_states(dftinfo.wave, emin, emax)
    super = Supercell(dftinfo.xtal, ORB_DICT)
    S = make_overlap_mat(occ_states)
    H = generate_H(super, DFTRAMO_EHT_PARAMS)
    
    if !isnothing(checkpoint)
        (psi_previous, num_electrons_left, num_raMO) = import_checkpoint(checkpoint)
    else
        num_electrons_left = sum([get(E_DICT, n.atom.name, 0) for n in dftinfo.xtal.atoms])
        num_raMO = 0
        psi_previous = diagm(ones(size(occ_states.coeff)[2]))
        psi_previous = ComplexF32.(repeat(psi_previous, 1, 1, 2)) #spin states to be implemented
    end
    
    low_psphere = Vector{Int}(undef, 0)
    next = iterate(runs)
    while next!==nothing
        (r, state) = next
        # print run information
        println(crayon"bold", "Run: ", crayon"light_cyan", r.name, crayon"!bold default")
        if r.type in keys(AO_RUNS)
            (low_psphere, psi_previous2, num_raMO2, num_electrons_left2) = loop_AO(
            super,
            r.sites,
            get(AO_RUNS, r.type, 0),
            num_electrons_left,
            num_raMO,
            r.name,
            DFTRAMO_EHT_PARAMS,
            occ_states,
            dftinfo.geo.basis,
            dftinfo.kpt,
            length.(collect.(dftinfo.wave.grange)),
            psi_previous,
            S,
            H,
            r.rsphere
            )
            site_list = r.sites # necessary for auto_psphere
        elseif r.type in CAGE_RUNS
            site_list = read_site_list(r.site_file)
            (low_psphere, psi_previous2, num_raMO2, num_electrons_left2) = loop_target_cluster_sp(
            super,
            site_list,
            r.radius,
            num_electrons_left,
            num_raMO,
            r.name,
            DFTRAMO_EHT_PARAMS,
            occ_states,
            dftinfo.geo.basis,
            dftinfo.kpt,
            length.(collect.(dftinfo.wave.grange)),
            psi_previous,
            S,
            H,
            r.rsphere
            )
        elseif r.type == "lcao"
            lcao_yaml = YAML.load_file(r.site_file)
            target = get(lcao_yaml, "target", nothing)
            isnothing(target) && error("Target cannot be blank for SALCs")
            for t in target
                !issubset(keys(t), keys(AO_RUNS)) && error("Target does not contain atomic orbitals. Please check.")
            end
            site_list = get(lcao_yaml, "lcao", nothing)
            if isa(site_list, Vector{Vector{Int}})
                for n in site_list
                    # TODO check the num atoms in lcao list match num targets
                    length(n) != length(target) && error("Mismatch between length of LCAO ", n, " and specified target.")
                end
            end
            (low_psphere, psi_previous2, num_raMO2, num_electrons_left2) = loop_LCAO(
            super,
            site_list,
            target,
            num_electrons_left,
            num_raMO,
            r.name,
            DFTRAMO_EHT_PARAMS,
            occ_states,
            dftinfo.geo.basis,
            dftinfo.kpt,
            length.(collect.(dftinfo.wave.grange)),
            psi_previous,
            S,
            H,
            r.rsphere
            )
        end
        # If auto_psphere is enabled, rerun and use the new run
        if auto_psphere && !isempty(low_psphere)
            # delete targets with low pspheres
            for i in reverse(low_psphere)
                deleteat!(site_list, i)
            end
            # set raMO analysis to where the first low psphere occurred
            # to prevent redundant raMO analysis
            deleteat!(site_list, collect(1:low_psphere[1]-1))
            e = num_electrons_left - (low_psphere[1]-1)*2
            raMO = num_raMO + (low_psphere[1]-1)
            (psi_previous, num_electrons_left, num_raMO) = import_checkpoint(string(r.name, "/", r.name, "_", raMO, "_", e, ".chkpt"))
            # If the last raMOs were the only ones with low_psphere, no need to rerun
            isempty(site_list) ? next = iterate(runs, state) : r.name = string(r.name, "_aps")
        else
            next = iterate(runs, state)
            num_electrons_left = num_electrons_left2
            num_raMO = num_raMO2
            psi_previous = psi_previous2
        end
    end
end

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
    warning("No energy units are specified. Defaulting to ", units, " (", software, ")")
    n = energy*e
    return n
end

"""
    parse_sites(sites::AbstractVector{<:AbstractString}) -> site_final::Vector{Int}

Parses a portion of the "sites" lines in the yaml file to return a Vector{Int} with valid
indices for targets.
e.g.,
```julia-repl
julia> parse_sites(["1:3", "3", "18:2:20"])
5-element Vector{Int64}:
  1
  2
  3
 18
 20
```
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

"""
    DFTraMO.parse_yaml_energy(x, i::InputOrigin)

Converts a string or number representing an energy input, possibly with units given, into an energy
in hartrees usable by DFTraMO. Units are not case-sensitive, and abbreviations are accepted as well
as full unit names.

DFT-raMO will automatically assume quantities without explicit unit specification match the units of
the software package which generated the input. If the software package is not specified or
supported, DFT-raMO will assume Hartree atomic units.
"""
parse_yaml_energy(x::Real, i::InputOrigin) = x * energy_conversion(i)

function parse_yaml_energy(x::AbstractString, i::InputOrigin)
    s = split(x)
    e = parse(Float64, first(s))
    # Strings which split into one substring likely don't contain units
    isone(length(s)) && return parse_yaml_energy(e, i)
    # Otherwise treat the units explicitly
    unit_str = s[2]
    conversion = if isempty(unit_str) || occursin(r"^au|ha|hartree"i, unit_str)
        1
    elseif occursin(r"^e(lectron-?)?v(olt)?"i, unit_str)
        Electrum.EV2HARTREE
    else
        error("Unit $unit_str was not recognized.")
    end
    return e * conversion
end

"""
    DFTraMO.parse_yaml_length(x, i::InputOrigin)

Converts a string or number representing an length input, possibly with units given, into a length
in bohr usable by DFTraMO. Units are not case-sensitive, and abbreviations are accepted as well as
full unit names.

Currently supported units are bohr (`bohr`, `a0`, or `au`), angstrom (`Å`, `ang`, or `angstrom`),
nanometer (`nm`), and picometer (`pm`).

DFT-raMO will automatically assume quantities without explicit unit specification match the units of
the software package which generated the input. If the software package is not specified or
supported, DFT-raMO will assume Hartree atomic units.
"""
parse_yaml_length(x::Real, i::InputOrigin) = x * length_conversion(i)

function parse_yaml_length(x::AbstractString, i::InputOrigin)
    s = split(x)
    e = parse(Float64, first(s))
    # Strings which split into one substring likely don't contain units
    isone(length(s)) && return parse_yaml_length(e, i)
    # Otherwise treat the units explicitly
    unit_str = s[2]
    conversion = if isempty(unit_str) || occursin(r"^(bohr|a[0u])"i, unit_str)
        1
    elseif occursin(r"^(Å|[AÅ]ng|[AÅ]ngstrom)"i, unit_str)
        Electrum.ANG2BOHR
    elseif unitstr == "nm"
        Electrum.ANG2BOHR * 10
    elseif unitstr == "pm"
        Electurm.ANG2BOHR / 100
    else
        error("Unit $unit_str was not recognized.")
    end
    return e * conversion
end

"""
    DFTraMO.read_yaml(file)

Reimplementation of `DFTraMO.read_run_yaml`.
"""
function read_yaml(io::IO)
    yaml = YAML.load(io)
    path = haskey(yaml, "path") ? yaml["path"] : "."
    origin = InputOrigin{Symbol(yaml["software"])}()
    dftinfo = raMOInput(path, origin)
    checkpoint = get(yaml, "checkpoint", nothing)
    # Check for empty string for robustness
    if isnothing(checkpoint) || isempty(checkpoint)
        @info "No checkpoint file found, starting new run from the beginning"
        checkpoint = ""
    elseif !isfile(checkpoint)
        error("Specified checkpoint $checkpoint is not a file.")
    else
        @info "Restarting calculation using checkpoint file $checkpoint"
    end
    # Use automatic conversion/type assertion for error
    auto_psphere::Bool = get(yaml, "auto_sphere", false)
    @info "Auto Psphere is " * (auto_psphere ? "en" : "dis") * "abled"
    # Get energy ranges (will need to perform unit inference)
    emin = get(yaml, "emin", nothing)
    emax = get(yaml, "emax", nothing)
    if isnothing(emin)
        emin = minimum(dftinfo.wave.energies)
        @info "No minimum energy specified: all bands below the maximum energy will be included."
    else
        emin = parse_yaml_energy(emin, origin)
    end
    if isnothing(emax)
        emax = fermi(dftinfo)
        @info "No maximum energy specified: all bands up to the Fermi energy will be included."
    else
        emax = parse_yaml_energy(emax, origin)
    end
    @assert emin < emax "emin ($emin Ha) is not less than emax ($emax Ha)!"
    @info string(
        "Energy range for reconstruction:\n",
        @sprintf("Minimum energy: %.3f Ha (%.3f eV)\n", emin, emin * Electrum.HARTREE2EV),
        @sprintf("Maximum energy: %.3f Ha (%.3f eV)\n", emax, emax * Electrum.HARTREE2EV)
    )
    # Get runs
    runs = get(data, "runs", nothing)
    @info "Performing $(length(runs)) runs."
    runlist = parse_runs(runs, dftinfo)
    return raMORuns(dftinfo, runlist, emin, emax, checkpoint, auto_psphere)
end

read_yaml(file) = open(read_yaml, file)
