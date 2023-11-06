"""
loop_target_cluster_sp()

Loops and runs DFT-raMO/Psphere analysis on a set of sp-based targets.
"""
function loop_target_cluster_sp(ramostatus::raMOStatus, sites)
    n = ramostatus.num_run
    ramoinput = ramostatus.ramoinput
    run = ramoinput[n]
    super = Supercell(ramoinput, ORB_DICT)
    # check to see which directory we're in
    if split(pwd(),'/')[end] != run.name && !isdir(run.name)
        mkdir(run.name)
    end
    cd(run.name) do
        if !isempty(readdir())
            @warn "Directory is not empty. Files will be deleted/overwritten."
            for n in readdir()
                rm(n, force=true)
            end
        end
        psphere = Vector{Float64}(undef, size(sites))
        remainders = Array{ComplexF32,1}(undef, 0)
        iter = ProgressBar(1:length(sites), unit="raMOs")
        for i in iter
            # check to see if we have electrons left
            if size(ramostatus.psi_previous)[2] == 0
                break
            end
            target = make_target_cluster_sp(sites, run.radius, i, super)
            (psi_previous2, psi_up, e_up, num_electrons_left2) = reconstruct_targets_DFT(target, DFTRAMO_EHT_PARAMS, ramostatus)
            isosurf = raMO_to_density(ramostatus.occ_states, psi_up, kptmesh(raMODFTData(ramoinput)), length.(collect.(PlanewaveWavefunction(ramoinput).grange)))
            (sphere, total, psphere[i]) = Psphere(RealDataGrid(real(isosurf), basis(super)), sites[i], run.rsphere)
            open(string(run.name, "_psphere_", run.rsphere, ".txt"), "a") do io
                print_psphere_terminal(iter, ramostatus.num_raMO+i, psphere[i], sites[i], io)
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
            remainders = ramostatus.psi_previous
        end
        cd("..")
        p = psphere_graph(psphere, ramostatus.num_raMO, run.rsphere); display(p)
        low_psphere = psphere_eval(psphere, super, sites)
        return (low_psphere, remainders, ramostatus.num_raMO+length(sites), ramostatus.num_electrons_left)
    end
end

function loop_AO(ramostatus::raMOStatus)
    n = ramostatus.num_run
    ramoinput = ramostatus.ramoinput
    run = ramoinput[n]
    super = Supercell(ramoinput, ORB_DICT)
    # check to see which directory we're in
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
            # check to see if we have electrons left
            if size(ramostatus.psi_previous)[2] == 0
                break
            end
            target = make_target_AO(run.sites[i], get(AO_RUNS, run.type, 0), super)
            (psi_previous2, psi_up, e_up, num_electrons_left2) = reconstruct_targets_DFT(target, DFTRAMO_EHT_PARAMS, ramostatus)
            isosurf = raMO_to_density(ramostatus.occ_states, psi_up, kptmesh(raMODFTData(ramoinput)), length.(collect.(PlanewaveWavefunction(ramoinput).grange)))
            pos = Vector(basis(super) * Electrum.BOHR2ANG * super[run.sites[i]].pos)
            (sphere, total, psphere[i]) = Psphere(RealDataGrid(real(isosurf), basis(super)), pos, run.rsphere)
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
            remainders = ramostatus.psi_previous
        end
        cd("..")
        p = psphere_graph(psphere, ramostatus.num_raMO, run.rsphere); display(p)
        # convert sites into cartesian; filter out atomic sites
        sites = [basis(super)*Electrum.BOHR2ANG*i.pos for i in PeriodicAtomList(super)][run.sites]
        low_psphere = psphere_eval(psphere, sites)
        return (low_psphere, remainders, ramostatus.num_raMO+length(run.sites), ramostatus.num_electrons_left)
    end
end

"""
loop_LCAO()

Loops and runs DFT-raMO/Psphere analysis on a set of LCAO targets.
"""
function loop_LCAO(ramostatus::raMOStatus)
    n = ramostatus.num_run
    ramoinput = ramostatus.ramoinput
    run = ramoinput[n]
    super = Supercell(ramoinput, ORB_DICT)
    lcao_yaml = YAML.load_file(run.site_file)
    target_orbital = get(lcao_yaml, "target", nothing)
    isnothing(target_orbital) && error("Target cannot be blank for SALCs")
    for t in target_orbital
        !issubset(keys(t), keys(AO_RUNS)) && error("Target does not contain atomic orbitals. Please check.")
    end
    site_list = get(lcao_yaml, "lcao", nothing)
    if isa(site_list, Vector{Vector{Int}})
        for n in site_list
            # TODO check the num atoms in lcao list match num targets
            length(n) != length(target_orbital) && error("Mismatch between length of LCAO ", n, " and specified target.")
        end
    end
    # check to see which directory we're in
    if split(pwd(),'/')[end] != run.name && !isdir(run.name)
        mkdir(run.name)
    end
    cd(run.name) do 
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
            # check to see if we have electrons left
            if size(ramostatus.psi_previous)[2] == 0
                break
            end
            target = make_target_lcao(site_list[i], target_orbital, super)
            (psi_previous2, psi_up, e_up, num_electrons_left2) = reconstruct_targets_DFT(target, DFTRAMO_EHT_PARAMS, ramostatus)
            isosurf = raMO_to_density(ramostatus.occ_states, psi_up, kptmesh(raMODFTData(ramoinput)), length.(collect.(PlanewaveWavefunction(ramoinput).grange)))
            pos = basis(super)*mp_lcao(site_list[i], PeriodicAtomList(super))*Electrum.BOHR2ANG
            (sphere, total, psphere[i]) = Psphere(RealDataGrid(real(isosurf), basis(super)), pos, run.rsphere)
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
            remainders = ramostatus.psi_previous
        end
        p = psphere_graph(psphere, ramostatus.num_raMO, run.rsphere); display(p)
        low_psphere = psphere_eval(psphere, super, site_list)
        return (low_psphere, remainders, ramostatus.num_raMO+length(site_list), ramostatus.num_electrons_left)
    end
end
