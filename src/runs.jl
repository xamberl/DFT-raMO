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

"""
print_psphere_terminal(iter, psphere, site)

Prints Psphere values to the terminal while raMO loops. The color of Psphere
    is arbitrarily color-coded.
    Green:        Psphere ≥ 0.5
    Yellow: 0.1 ≤ Psphere < 0.5
    Red:          Psphere < 0.1
"""
function print_psphere_terminal(iter, num_raMO, psphere, site)
    col = Crayon(foreground = :green)
    psphere < 0.5 ? col = Crayon(foreground = :light_yellow) : nothing
    psphere < 0.1 ? col = Crayon(foreground = :light_red) : nothing
    println(
        iter,
        num_raMO,
        ", Psphere: ",
        col, @sprintf("%.3f", psphere),
        Crayon(foreground = :default),
        @sprintf(" at site [%.3f, %.3f, %.3f]", site[1], site[2], site[3]))
end
    
function psphere_graph(psphere::Vector{Float64}, num_raMO::Int, rsphere::Float64)
    p = lineplot(
        collect(1:length(psphere)).+num_raMO,
        psphere,
        title="Psphere",
        name=string("rsphere @ ", rsphere),
        xlabel="raMO",
        ylabel="Psphere",
        ylim=(0,1))
    return p
end