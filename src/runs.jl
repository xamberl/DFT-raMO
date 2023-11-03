"""
    loop_target_cluster_sp()

Loops and runs DFT-raMO/Psphere analysis on a set of sp-based targets.
"""
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
    rsphere,
    discard
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
        remainders = Array{ComplexF32,1}(undef, 0)
        iter = ProgressBar(1:length(voids_list), unit="raMOs")
        for i in iter
            target = make_target_cluster_sp(voids_list, void_radius, i, super)
            (psi_previous2, psi_up, e_up, num_electrons_left2) = reconstruct_targets_DFT(target, ehtparams, ramostatus)
            isosurf = raMO_to_density(occ_states, psi_up, kpt, grange)
            (sphere, total, psphere[i]) = Psphere(RealDataGrid(real(isosurf),super.atomlist.basis), voids_list[i], rsphere)
            open(string(run_name, "_psphere_", rsphere, ".txt"), "a") do io
                print_psphere_terminal(iter, num_raMO+i, psphere[i], voids_list[i], io)
            end
            # Check if we are in discard ramo mode
            if discard
                # print out raMO but not checkpoint
                write_to_XSF(isosurf, super.atomlist, string(run_name, "_", i, ".xsf"))
                open(string(run_name, "_", i, ".raMO"), "w") do io
                    write(io, psi_up)
                end
            else
                # update remainders and number of electrons left
                psi_previous = psi_previous2; num_electrons_left = Int(num_electrons_left2)
                output_files(run_name, num_electrons_left, num_raMO+i, super, isosurf, psi_previous, psi_up)
            end
            remainders = psi_previous
        end
        cd("..")
        p = psphere_graph(psphere, num_raMO, rsphere); display(p)
        low_psphere = psphere_eval(psphere, super, voids_list)
        return (low_psphere, remainders, num_raMO+length(voids_list), num_electrons_left)
    end
end

#=="""
    loop_AO()

Loops and runs DFT-raMO/Psphere analysis on a set of atomic orbital targets.
"""
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
    rsphere,
    discard
) ## MUTABLE STRUCT RECOMMENDED
    # check to see which directory we're in
    if split(pwd(),'/')[end] != run_name && !isdir(run_name)
        mkdir(run_name)
    end
    cd(run_name) do 
        # check to see if directory is empty. if not, send warning before deleting
        # TODO: in case someone puts a file in there, just wipe the necessary files rather
        # than wiping everything.
        if !isempty(readdir())
            @warn "Directory is not empty. Files will be deleted/overwritten."
            for n in readdir()
                rm(n, force=true)
            end
        end
        psphere = Vector{Float64}(undef, size(atom_list))
        remainders = Array{ComplexF32,1}(undef, 0) # TODO
        iter = ProgressBar(eachindex(atom_list), unit="raMOs")
        for i in iter
            target = make_target_AO(atom_list[i], target_orbital, super)
            (psi_previous2, psi_up, e_up, num_electrons_left2) = reconstruct_targets_DFT(
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
            isosurf = raMO_to_density(occ_states, psi_up, kpt, grange)
            pos = Vector(basis(super) * Electrum.BOHR2ANG * super[atom_list[i]].pos)
            (sphere, total, psphere[i]) = Psphere(RealDataGrid(real(isosurf),basis(super)), pos, rsphere)
            open(string(run_name, "_psphere_", rsphere, ".txt"), "a") do io
                print_psphere_terminal(iter, num_raMO+i, psphere[i], pos, io)
            end
            # Check if we are in discard ramo mode
            if discard
                # print out raMO but not checkpoint
                write_to_XSF(isosurf, super.atomlist, string(run_name, "_", i, ".xsf"))
                open(string(run_name, "_", i, ".raMO"), "w") do io
                    write(io, psi_up)
                end
            else
                # update remainders and number of electrons left
                psi_previous = psi_previous2; num_electrons_left = Int(num_electrons_left2)
                output_files(run_name, num_electrons_left, num_raMO+i, super, isosurf, psi_previous, psi_up)
            end
            remainders = psi_previous
        end
        cd("..")
        p = psphere_graph(psphere, num_raMO, rsphere); display(p)
        low_psphere = psphere_eval(psphere, super, atom_list)
        return (low_psphere, remainders, num_raMO+length(atom_list), num_electrons_left)
    end
end==#

function loop_AO(ramostatus::raMOStatus)
    # check to see which directory we're in
    n = ramostatus.num_run
    ramoinput = ramostatus.ramoinput
    run = ramoinput[n]
    super = Supercell(ramoinput, ORB_DICT)
    if split(pwd(),'/')[end] != run.name && !isdir(run.name)
        mkdir(run.name)
    end
    cd(run.name) do 
        # check to see if directory is empty. if not, send warning before deleting
        # TODO: in case someone puts a file in there, just wipe the necessary files rather
        # than wiping everything.
        if !isempty(readdir())
            @warn "Directory is not empty. Files will be deleted/overwritten."
            for n in readdir()
                rm(n, force=true)
            end
        end
        psphere = Vector{Float64}(undef, size(run.sites))
        remainders = Array{ComplexF32,1}(undef, 0) # TODO
        iter = ProgressBar(eachindex(run.sites), unit="raMOs")
        for i in iter
            target = make_target_AO(run.sites[i], get(AO_RUNS, run.type, 0), super)
            (psi_previous2, psi_up, e_up, num_electrons_left2) = reconstruct_targets_DFT(target, DFTRAMO_EHT_PARAMS, ramostatus)
            isosurf = raMO_to_density(ramostatus.occ_states, psi_up, kptmesh(raMODFTData(ramoinput)), length.(collect.(PlanewaveWavefunction(ramoinput).grange)))
            pos = Vector(basis(super) * Electrum.BOHR2ANG * super[run.sites[i]].pos)
            (sphere, total, psphere[i]) = Psphere(RealDataGrid(real(isosurf),basis(super)), pos, run.rsphere)
            open(string(run.name, "_psphere_", run.rsphere, ".txt"), "a") do io
                print_psphere_terminal(iter, ramostatus.num_raMO+i, psphere[i], pos, io)
            end
            # Check if we are in discard ramo mode
            if ramoinput.discard
                # print out raMO but not checkpoint
                write_to_XSF(isosurf, super.atomlist, string(run.name, "_", i, ".xsf"))
                open(string(run.name, "_", i, ".raMO"), "w") do io
                    write(io, psi_up)
                end
            else
                # update remainders and number of electrons left
                ramostatus.psi_previous = psi_previous2
                ramostatus.num_electrons_left = Int(num_electrons_left2)
                output_files(run.name, ramostatus.num_electrons_left, ramostatus.num_raMO+i, super, isosurf, ramostatus.psi_previous, psi_up)
            end
        end
        cd("..")
        p = psphere_graph(psphere, ramostatus.num_raMO, run.rsphere); display(p)
        low_psphere = psphere_eval(psphere, super, run.sites)
        return (low_psphere, remainders, ramostatus.num_raMO+length(run.sites), num_electrons_left)
    end
end

"""
    loop_LCAO()

Loops and runs DFT-raMO/Psphere analysis on a set of LCAO targets.
"""
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
    rsphere,
    discard
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
        remainders = Array{ComplexF32,1}(undef, 0)
        iter = ProgressBar(1:length(site_list), unit="raMOs")
        for i in iter
            target = make_target_lcao(site_list[i], target_orbital, super)
            (psi_previous2, psi_up, e_up, num_electrons_left2) = reconstruct_targets_DFT(target, ehtparams, ramostatus)
            isosurf = raMO_to_density(occ_states, psi_up, kpt, grange)
            pos = super.atomlist.basis*mp_lcao(site_list[i], super.atomlist)*Electrum.BOHR2ANG
            (sphere, total, psphere[i]) = Psphere(RealDataGrid(real(isosurf),super.atomlist.basis), pos, rsphere)
            open(string(run_name, "_psphere_", rsphere, ".txt"), "a") do io
                print_psphere_terminal(iter, num_raMO+i, psphere[i], pos, io)
            end
            # Check if we are in discard ramo mode
            if discard
                # print out raMO but not checkpoint
                write_to_XSF(isosurf, super.atomlist, string(run_name, "_", i, ".xsf"))
                open(string(run_name, "_", i, ".raMO"), "w") do io
                    write(io, psi_up)
                end
            else
                # update remainders and number of electrons left
                psi_previous = psi_previous2; num_electrons_left = Int(num_electrons_left2)
                output_files(run_name, num_electrons_left, num_raMO+i, super, isosurf, psi_previous, psi_up)
            end
            remainders = psi_previous
        end
        cd("..")
        p = psphere_graph(psphere, num_raMO, rsphere); display(p)
        low_psphere = psphere_eval(psphere, super, site_list)
        return (low_psphere, remainders, num_raMO+length(site_list), num_electrons_left)
    end
end
