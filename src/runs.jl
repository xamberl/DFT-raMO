function loop_target_cluster_sp(
    super,
    voids_list,
    void_radius,
    num_electrons_left,
    run_name,
    ehtparams,
    occ_states,
    geo_basis,
    kpt,
    psi_previous,
    S,
    H,
    rsphere,
    )
    psphere = Vector{Float64}(undef, size(voids_list))
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
        isosurf = psi_to_isosurf(5000000, super.atomlist, kpt, occ_states, psi_up)
        (sphere, total, psphere[i]) = Psphere(RealDataGrid(real(isosurf),super.atomlist.basis), voids_list[i], rsphere)
        println(iter, string("Psphere: ", psphere[i], " at site ", voids_list[i]))
        write_to_XSF(isosurf, super.atomlist, string(run_name, "_", i, "_", num_electrons_left, ".xsf"))
        # for now, write only one spin as save state
        writedlm(string(run_name, "_psi_prev_", i, "_", num_electrons_left, ".txt"), psi_previous[:,:,1])
        # same with the raMO function itself
        writedlm(string(run_name, "_psi_", i, "_", num_electrons_left, ".txt"), psi_up)
    end
    writedlm(string(run_name, "_psphere_", rsphere, ".txt"), psphere)
    lineplot(collect(1:length(psphere)), psphere, title="Psphere", name=string("rsphere @ ", rsphere), xlabel="raMO", ylabel="Psphere", ylim=(0,1))
end