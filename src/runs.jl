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
    cd(run_name)
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
        isosurf = psi_to_isosurf2(occ_states, psi_up, kpt, grange)
        (sphere, total, psphere[i]) = Psphere(RealDataGrid(real(isosurf),super.atomlist.basis), voids_list[i], rsphere)
        print_psphere_terminal(iter, num_raMO+i, psphere[i], voids_list[i])
        output_files(run_name, num_electrons_left, num_raMO+i, super, isosurf, psi_previous, psi_up)
        remainders = psi_previous
    end
    writedlm(string(run_name, "_psphere_", rsphere, ".txt"), psphere)
    cd("..")
    p = psphere_graph(psphere, num_raMO, rsphere); display(p)
    low_psphere = psphere_eval(psphere, super, voids_list)
    return (low_psphere, remainders, num_raMO+length(voids_list), num_electrons_left)
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
    cd(run_name)
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
        isosurf = psi_to_isosurf2(occ_states, psi_up, kpt, grange)
        pos = Vector(super.atomlist.basis*Electrum.BOHR2ANG*super.atomlist[atom_list[i]].pos)
        (sphere, total, psphere[i]) = Psphere(RealDataGrid(real(isosurf),super.atomlist.basis), pos, rsphere)
        print_psphere_terminal(iter, num_raMO+i, psphere[i], pos)
        output_files(run_name, num_electrons_left, num_raMO+i, super, isosurf, psi_previous, psi_up)
        remainders = psi_previous
    end
    writedlm(string(run_name, "_psphere_", rsphere, ".txt"), psphere)
    cd("..")
    p = psphere_graph(psphere, num_raMO, rsphere); display(p)
    low_psphere = psphere_eval(psphere, super, atom_list)
    return (low_psphere, remainders, num_raMO+length(atom_list), num_electrons_left)
end

function dftramo_run(filename::AbstractString)
    (runs, checkpoint, auto_psphere, dftinfo) = read_run_yaml(filename)
    ehtparams = read_eht_params("DFT_raMO_eht_parms.dat")
    occ_states = get_occupied_states(dftinfo.wave, fermi.fermi*Electrum.EV2HARTREE)
    super = Supercell(dftramo.xtal, DFTraMO.orb_dict)
    S = make_overlap_mat(occ_states)
    H = DFTraMO.generate_H(super, ehtparams)

    if !isnothing(checkpoint)
        (psi_previous, num_electrons_left, num_raMO) = import_psi_previous(checkpoint)
    else
        num_electrons_left = sum([get(e_dict, n.atom.name, 0) for n in dftinfo.xtal.atoms])
        num_raMO = 0
        psi_previous = diagm(ones(size(occ_states.coeff)[2]))
        psi_previous = ComplexF32.(repeat(psi_previous, 1, 1, 2))
    end

    # This loop needs some design... WIP
    low_psphere = Vector{Int}(undef, 0)
    auto_psphere ? aps = "_aps" : aps = ""
    next = iterate(runs)
    while next!==nothing
    (r, state) = next
        if in(r.type, AO_RUNS)
            (low_psphere, psi_previous2, num_raMO2, num_electrons_left2) = DFTraMO.loop_AO(
                super,
                r.sites,
                get(AO_RUNS, r.type, 0),
                num_electrons_left,
                num_raMO,
                string(r.name, aps),
                ehtparams,
                occ_states,
                dftinfo.geo.basis,
                dftinfo.kpt,
                length.(collect.(dftinfo.wave.grange)),
                psi_previous,
                S,
                H,
                r.rsphere
                )
        elseif in(r.type, CAGE_RUNS)
            site_list = read_site_list(r.site_file)
            (low_psphere, psi_previous2, num_raMO2, num_electrons_left2) = loop_target_cluster_sp(
                super,
                site_list,
                r.radius,
                num_electrons_left,
                num_raMO,
                string(r.name, aps),
                ehtparams,
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
        if auto_psphere
           for i in reverse(low_psphere)
                deleteat!(atom_list, i)
           end
        else
            next = iterate(runs, state)
        end
    end
end