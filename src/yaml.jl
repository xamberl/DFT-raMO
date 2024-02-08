"""
    dftramo_run(filename::AbstractString)

Automatically runs DFTraMO from a configuration yaml file.
"""
function dftramo_run(filename::AbstractString)
    # read inputs
    ramoinput = read_yaml(filename)
    # initialize data
    ramostatus = raMOStatus(ramoinput)
    
    next = iterate(ramoinput)
    while next!==nothing
        (r, state) = next
        ramostatus.num_run = state-1
        # check to see if we have states left
        size_basis(ramostatus) == 0 && info("No more states in the basis set. Stopping DFT-raMO.") && break
        # print run information
        @info "Run: $(r.name)"

        # go into directory
        split(wd(ramoinput),'/')[end] != r.name && !isdir(r.name) && mkdir(r.name)
        cd(r.name) do
            if !isempty(readdir())
                @warn "Directory is not empty. Files will be deleted/overwritten."
                [rm(n, force=true, recursive=true) for n in readdir()]
            end
            psphere = Vector{Float64}(undef, length(r.sites))
            psphere_sites = Vector{SVector{3, Float64}}(undef, length(r.sites))
            ramostatus = run_ramo(psphere, psphere_sites, r, ramostatus)
            # If auto_psphere is enabled, rerun and use the new run
            low_psphere = psphere_eval(psphere, psphere_sites)
            if mode(ramoinput) == "auto_psphere" && length(low_psphere)!=0
                # delete targets with low pspheres
                for i in reverse(low_psphere)
                    deleteat!(psphere_sites, i)
                end
                # set raMO analysis to where the first low psphere occurred
                # to prevent redundant raMO analysis
                deleteat!(psphere_sites, collect(1:low_psphere[1]-1))
                raMO = length(psphere) - low_psphere[1]
                e = ramostatus.num_electrons_left + raMO*2
                (ramostatus.psi_previous, ramostatus.num_electrons_left, ramostatus.num_raMO) = import_checkpoint(string(r.name, "_", ramostatus.num_raMO-raMO-1, "_", e+2, ".chkpt"))
                # If the last raMOs were the only ones with low_psphere, no need to rerun
                if isempty(psphere_sites)
                    next = iterate(ramoinput, state)
                else
                    mkdir("aps")
                    cd("aps") do 
                        psphere = Vector{Float64}(undef, length(psphere_sites))
                        psphere_sites = Vector{SVector{3, Float64}}(undef, length(psphere_sites))
                        ramostatus = run_ramo(psphere, psphere_sites, r, ramostatus, low_psphere=low_psphere)
                    end
                    next = iterate(ramoinput, state)
                end
            else
                next = iterate(ramoinput, state)
            end
        end
    end
end

"""
    run_ramo(psphere, psphere_sites, r, ramostatus; low_psphere = Vector{Int}(undef, 0))

Runs the basic raMO code (necessary in a separate function for auto_psphere functionality)
"""
function run_ramo(psphere, psphere_sites, r, ramostatus; low_psphere = Vector{Int}(undef, 0))
    iter = ProgressBar(1:length(psphere), unit="raMOs")
    for i in iter
        # check to see if we have electrons left
        if size(ramostatus.psi_previous)[2] == 0
            # truncate psphere and psphere_sites
            if i == 1
                psphere = Vector{Float64}(undef, 0)
                psphere_sites = Vector{SVector{3, Float64}}(undef, 0)
            else
                psphere = psphere[1:i-1]
                psphere_sites = psphere_sites[1:i-1]
            end
            break
        end
        displacement = [0,0,0]
        # generate target
        if r.type in keys(AO_RUNS)
            target_sites = copy(r.sites)
            for n in reverse(low_psphere)
                deleteat!(target_sites, n)
            end
            length(low_psphere) > 0 && deleteat!(target_sites, collect(1:low_psphere[1]-1)) 
            target = make_target_AO(target_sites[i], get(AO_RUNS, r.type, 0), ramostatus.supercell)
            if norm(r.direction) > 0
                # displacement in fractional coordinates
                displacement = inv(basis(raMODFTData(raMOInput(ramostatus))))*(r.direction * r.radius * Electrum.ANG2BOHR)
            end
            sites = basis(ramostatus.supercell)*Electrum.BOHR2ANG*PeriodicAtomList(ramostatus.supercell)[target_sites[i]].pos+(r.direction * r.radius)
        elseif r.type in CAGE_RUNS # for now we only have one type
            sites = read_site_list(joinpath(wd(raMOInput(ramostatus)), r.site_file))
            for n in reverse(low_psphere)
                deleteat!(sites, n)
            end
            length(low_psphere) > 0 && deleteat!(sites, collect(1:low_psphere[1]-1)) 
            target = make_target_cluster_sp(sites, r.radius, i, ramostatus.supercell)
            sites = sites[i]
        elseif r.type == "lcao"
            lcao_yaml = YAML.load_file(joinpath(wd(raMOInput(ramostatus)), r.site_file))
            (target_orbital, site_list) = target_lcao(lcao_yaml)
            for n in reverse(low_psphere)
                deleteat!(site_list, n)
            end
            length(low_psphere) > 0 && deleteat!(site_list, collect(1:low_psphere[1]-1)) 
            target = make_target_lcao(site_list[i], target_orbital, ramostatus.supercell)
            sites = basis(ramostatus.supercell)*mp_lcao(site_list[i], PeriodicAtomList(ramostatus.supercell))*Electrum.BOHR2ANG
        end
        (psi_previous2, psi_up, e_up, num_electrons_left2) = reconstruct_targets_DFT(target, DFTRAMO_EHT_PARAMS, ramostatus, displacement)
        isosurf = raMO_to_density(OccupiedStates(ramostatus), psi_up, kptmesh(raMODFTData(ramostatus)), length.(collect.(PlanewaveWavefunction(raMOInput(ramostatus)).grange)))
        (sphere, total, psphere[i]) = Psphere(RealDataGrid(real(isosurf), basis(ramostatus.supercell)), sites, r.rsphere)
        psphere_sites[i] = sites
        open(string(r.name, "_psphere_", r.rsphere, ".txt"), "a") do io
            print_psphere_terminal(iter, ramostatus.num_raMO+1, psphere[i], sites, io)
        end
        # Check if we are in discard ramo mode
        if mode(raMOInput(ramostatus)) == "discard"
            # print out raMO but not checkpoint
            write_to_XSF(isosurf, PeriodicAtomList(ramostatus.supercell), string(r.name, "_", i, ".xsf"))
            open(string(r.name, "_", i, ".raMO"), "w") do io
                write(io, psi_up)
            end
        else
            # update remainders and number of electrons left
            ramostatus.psi_previous = psi_previous2
            ramostatus.num_electrons_left = Int(num_electrons_left2)
            ramostatus.num_raMO += 1
            output_files(r.name, ramostatus.num_electrons_left, ramostatus.num_raMO, ramostatus.supercell, isosurf, ramostatus.psi_previous, psi_up)
        end
    end
    return ramostatus
end

"""
    target_lcao(yaml::Dict) -> target_orbital

Parses the yaml file for LCAOs.
"""
function target_lcao(yaml::Dict)
    target_orbital = get(yaml, "target", nothing)
    isnothing(target_orbital) && error("Target cannot be blank for LCAOs")
    for t in target_orbital
        !issubset(keys(t), keys(DFTraMO.AO_RUNS)) && error("Target does not contain atomic orbitals. Please check.")
    end
    site_list = get(yaml, "lcao", nothing)
    if isa(site_list, Vector{Vector{Int}})
        for n in site_list
            # TODO check the num atoms in lcao list match num targets
            length(n) != length(target_orbital) && error("Mismatch between length of LCAO ", n, " and specified target.")
        end
    end
    return (target_orbital, site_list)
end

"""
    parse_runs(runs::Vector{Dict{Any, Any}}, dftinfo::raMODFTData, origin::InputOrigin) -> runlist::Vector{RunInfo}

Parse the run section of the input yaml file and returns Vector{RunInfo}.
"""
function parse_runs(runs::Vector{Dict{Any, Any}}, dftinfo::raMODFTData, origin::InputOrigin)
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
        # check for displaced tag
        if split(type)[1] == "displaced"
            displaced = true
            type = split(type)[2]
        else
            displaced = false
        end
        !in(type, union(keys(AO_RUNS), CAGE_RUNS, ["lcao"])) && error("type ", type, " is an invalid entry.")
        println("   type: ", cr_b, type, cr_d)

        site_file = get(runs[n], "site_file", "")

        # site_file is not needed for an AO run
        if in(type, keys(AO_RUNS))
            println("   AO type run. site_file will be ignored.")
            site_file = ""
            site_list = Electrum.supercell(dftinfo)
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
            if in(type, keys(AO_RUNS)) && !isdigit(sites[1][1])
                e = string(sites[1])
                !in(e, Electrum.ELEMENTS) && error(e, " is not a valid element.")
                valid_atoms = findall(x -> Electrum.name(x) == e, Electrum.supercell(dftinfo))
                isempty(valid_atoms) && error("Element ", e, "was not found in the system.") 
                length(sites) > 1 ? sites_final = valid_atoms[parse_sites(sites[2:end])] : sites_final = valid_atoms
            else
                sites_final = parse_sites(sites)
                length(sites_final) > length(site_list) && error("Number of sites specified exceeds site_list length.")
            end
        end
        println("   sites: ", cr_b, get(runs[n], "sites", nothing), cr_d)

        # get direction for displaced AOs
        if displaced
            direction = get(runs[n], "direction", [0,0,0])
            direction == [0,0,0] && @error "Please specify a nonzero vector."
            !(isa(direction, Vector{<:Number}) && length(direction) == 3) && @error "direction is not a valid 3-element vector."
            direction = direction/norm(direction) # normalize direction
        else
            cart = true
            direction = [0,0,0]
        end
        
        radius = get(runs[n], "radius", 0.0)
        if isnothing(radius)
            radius = 0.0
        end
        if in(type, keys(AO_RUNS))
            if displaced
                radius = parse_yaml_length(radius, origin)
                iszero(radius) && @error "Please specify a nonzero radius."
                println("   radius: ", cr_b, radius, " Å", cr_d)
            else
                println("   Radius is ignored for atomic orbital type runs.")
            end
        else
            radius = parse_yaml_length(radius, origin)
            println("   radius: ", cr_b, radius, " Å", cr_d)
        end

        rsphere = get(runs[n], "rsphere") do
            println(cr_y, "   rsphere was not specified. Default value of 3.0 Å is applied.", cr_d)
            return 3.0
        end
        rsphere = parse_yaml_length(rsphere, origin)
        println("   rsphere: ", cr_b, rsphere, " Å", cr_d)

        runlist[n] = RunInfo(name, type, site_file, sites_final, direction, radius, rsphere)
    end
    return runlist
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
    conversion = if isempty(unit_str) || occursin(r"^(Å|[AÅ]ng|[AÅ]ngstrom)"i, unit_str)
        1
    elseif occursin(r"^(bohr|a[0u])"i, unit_str)
        Electrum.BOHR2ANG
    elseif unitstr == "nm"
        10
    elseif unitstr == "pm"
        0.01
    else
        error("Unit $unit_str was not recognized.")
    end
    return e * conversion
end

"""
    DFTraMO.read_yaml(file)

Reads an input yaml and returns an raMOInput object.
"""
function read_yaml(io::IO)
    yaml = YAML.load(io)
    path = haskey(yaml, "path") ? yaml["path"] : "."
    origin = InputOrigin{Symbol(yaml["software"])}()
    dftinfo = raMODFTData(path, origin)
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
    # check for mode
    mode = get(yaml, "mode", "default")
    in(lowercase(mode), ["default", "discard", "auto_psphere"]) ? mode = lowercase(mode) : error("mode must be 'default', 'discard', or 'auto_psphere'.")
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
    runs = get(yaml, "runs", nothing)
    @info "Performing $(length(runs)) runs."
    runlist = parse_runs(runs, dftinfo, origin)
    return raMOInput(dftinfo, runlist, emin, emax, checkpoint, mode, pwd())
end

read_yaml(file) = open(read_yaml, file)
