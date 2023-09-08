function loop_target_cluster_sp(
    super,
    voids_list,
    void_radius,
    num_electrons_left,
    num_raMO,
    run_name,
    ehtparams,
    occ_states,
    geo_basis,
    kpt,
    grange,
    psi_previous,
    S,
    H,
    rsphere
)
    # check to see which directory we're in
    if split(pwd(),'/')[end] != run_name && !isdir(run_name)
        mkdir(run_name)
    end
    cd(run_name) do
        if !isempty(readdir())
            @warn "Directory is not empty. Files will be deleted/overwritten."
            for n in readdir()
                rm(n, force=true)
            end
        end
        psphere = Vector{Float64}(undef, size(voids_list))
        remainders = []
        iter = ProgressBar(1:length(voids_list), unit="raMOs")
        for i in iter
            target = make_target_cluster_sp(voids_list, void_radius, i, super)
            (psi_previous, psi_up, e_up, num_electrons_left) = reconstruct_targets_DFT(
                target,
                num_electrons_left,
                run_name,
                super,
                ehtparams,
                occ_states,
                geo_basis,
                kpt,
                psi_previous,
                S,
                H,
                false,
                "",
            )
            isosurf = psi_to_isosurf(occ_states, psi_up, kpt, grange)
            (sphere, total, psphere[i]) = Psphere(RealDataGrid(real(isosurf),super.atomlist.basis), voids_list[i], rsphere)
            open(string(run_name, "_psphere_", rsphere, ".txt"), "a") do io
                print_psphere_terminal(iter, num_raMO+i, psphere[i], voids_list[i], io)
            end
            output_files(run_name, num_electrons_left, num_raMO+i, super, isosurf, psi_previous, psi_up)
            remainders = psi_previous
        end
        cd("..")
        p = psphere_graph(psphere, num_raMO, rsphere); display(p)
        low_psphere = psphere_eval(psphere, super, voids_list)
        return (low_psphere, remainders, num_raMO+length(voids_list), num_electrons_left)
    end
end

function loop_AO(
    super,
    atom_list,
    target_orbital,
    num_electrons_left,
    num_raMO,
    run_name,
    ehtparams,
    occ_states,
    geo_basis,
    kpt,
    grange,
    psi_previous,
    S,
    H,
    rsphere
)
    # check to see which directory we're in
    if split(pwd(),'/')[end] != run_name && !isdir(run_name)
        mkdir(run_name)
    end
    cd(run_name) do 
        # check to see if directory is empty. if not, send warning before deleting
        if !isempty(readdir())
            @warn "Directory is not empty. Files will be deleted/overwritten."
            for n in readdir()
                rm(n, force=true)
            end
        end
        psphere = Vector{Float64}(undef, size(atom_list))
        remainders = []
        iter = ProgressBar(1:length(atom_list), unit="raMOs")
        for i in iter
            target = make_target_AO(atom_list[i], target_orbital, super)
            (psi_previous, psi_up, e_up, num_electrons_left) = reconstruct_targets_DFT(
                target,
                num_electrons_left,
                run_name,
                super,
                ehtparams,
                occ_states,
                geo_basis,
                kpt,
                psi_previous,
                S,
                H,
                false,
                "",
            )
            isosurf = psi_to_isosurf(occ_states, psi_up, kpt, grange)
            pos = Vector(super.atomlist.basis*Electrum.BOHR2ANG*super.atomlist[atom_list[i]].pos)
            (sphere, total, psphere[i]) = Psphere(RealDataGrid(real(isosurf),super.atomlist.basis), pos, rsphere)
            open(string(run_name, "_psphere_", rsphere, ".txt"), "a") do io
                print_psphere_terminal(iter, num_raMO+i, psphere[i], pos, io)
            end
            output_files(run_name, num_electrons_left, num_raMO+i, super, isosurf, psi_previous, psi_up)
            remainders = psi_previous
        end
        cd("..")
        p = psphere_graph(psphere, num_raMO, rsphere); display(p)
        low_psphere = psphere_eval(psphere, super, atom_list)
        return (low_psphere, remainders, num_raMO+length(atom_list), num_electrons_left)
    end
end

function loop_LCAO(
    super,
    site_list, #atom_list,
    target_orbital,
    num_electrons_left,
    num_raMO,
    run_name,
    ehtparams,
    occ_states,
    geo_basis,
    kpt,
    grange,
    psi_previous,
    S,
    H,
    rsphere
)
    # check to see which directory we're in
    if split(pwd(),'/')[end] != run_name && !isdir(run_name)
        mkdir(run_name)
    end
    cd(run_name) do 
        if !isempty(readdir())
            @warn "Directory is not empty. Files will be deleted/overwritten."
            for n in readdir()
                rm(n, force=true)
            end
        end
        psphere = Vector{Float64}(undef, length(site_list))
        remainders = []
        iter = ProgressBar(1:length(site_list), unit="raMOs")
        for i in iter
            target = make_target_lcao(site_list[i], target_orbital, super)
            (psi_previous, psi_up, e_up, num_electrons_left) = reconstruct_targets_DFT(
                target,
                num_electrons_left,
                run_name,
                super,
                ehtparams,
                occ_states,
                geo_basis,
                kpt,
                psi_previous,
                S,
                H,
                false,
                "",
            )
            isosurf = psi_to_isosurf(occ_states, psi_up, kpt, grange)
            pos = super.atomlist.basis*mp_lcao(site_list[i], super.atomlist)*Electrum.BOHR2ANG
            #pos = Vector(super.atomlist.basis*Electrum.BOHR2ANG*super.atomlist[site_list[i][1]].pos)
            (sphere, total, psphere[i]) = Psphere(RealDataGrid(real(isosurf),super.atomlist.basis), pos, rsphere)
            open(string(run_name, "_psphere_", rsphere, ".txt"), "a") do io
                print_psphere_terminal(iter, num_raMO+i, psphere[i], pos, io)
            end
            output_files(run_name, num_electrons_left, num_raMO+i, super, isosurf, psi_previous, psi_up)
            remainders = psi_previous
        end
        cd("..")
        p = psphere_graph(psphere, num_raMO, rsphere); display(p)
        low_psphere = psphere_eval(psphere, super, site_list)
        return (low_psphere, remainders, num_raMO+length(site_list), num_electrons_left)
    end
end

"""
    dftramo_run(filename::AbstractString)

Reads a run.yaml file and executes DFTraMO.
"""
function dftramo_run(filename::AbstractString, software::AbstractString="vasp")
    #redirect_stdio(stdout=string(splitext(filename)[1], ".log")) do
        (runs, checkpoint, auto_psphere, dftinfo) = read_run_yaml(filename, software)
        occ_states = get_occupied_states(dftinfo.wave, dftinfo.fermi.fermi*Electrum.EV2HARTREE)
        super = Supercell(dftinfo.xtal, orb_dict)
        S = make_overlap_mat(occ_states)
        H = DFTraMO.generate_H(super, DFTRAMO_EHT_PARAMS)
        
        if !isnothing(checkpoint)
            (psi_previous, num_electrons_left, num_raMO) = import_chkpt(checkpoint)
        else
            num_electrons_left = sum([get(e_dict, n.atom.name, 0) for n in dftinfo.xtal.atoms])
            num_raMO = 0
            psi_previous = diagm(ones(size(occ_states.coeff)[2]))
            psi_previous = ComplexF32.(repeat(psi_previous, 1, 1, 2))
        end
        
        low_psphere = Vector{Int}(undef, 0)
        next = iterate(runs)
        while next!==nothing
            (r, state) = next
            # print run information
            println(crayon"bold", "Run: ", crayon"light_cyan", r.name, crayon"!bold default")
            if r.type in keys(AO_RUNS)
                (low_psphere, psi_previous2, num_raMO2, num_electrons_left2) = DFTraMO.loop_AO(
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
                (psi_previous, num_electrons_left, num_raMO) = import_chkpt(string(r.name, "/", r.name, "_", raMO, "_", e, ".chkpt"))
                # If the last raMOs were the only ones with low_psphere, no need to rerun
                isempty(site_list) ? next = iterate(runs, state) : r.name = string(r.name, "_aps")
            else
                next = iterate(runs, state)
                num_electrons_left = num_electrons_left2
                num_raMO = num_raMO2
                psi_previous = psi_previous2
            end
        end
    #end
end
