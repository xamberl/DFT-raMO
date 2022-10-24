"""
get_eht_params(atom_num::Int, eht_params::ehtParams) -> OrbitalParams

Search for parameters in the loaded ehtParams for corresponding atom.
"""
function get_eht_params(atom_num::Int, eht_params::ehtParams)
    return eht_params.data[atom_num, :]
end


"""
    make_overlap_mat(occ_coeff, kpoint_repeating::Vector{Bool})

    Creates overlap matrix with occupied coefficients
    If the kpoints are the same, they can overlap. Otherwise, they do not.
"""
function make_overlap_mat(occ_coeff, kpoint_repeating::Vector{Bool})
    num_occ_states = length(kpoint_repeating)
    S = zeros(ComplexF64,num_occ_states,num_occ_states)
    # Creates list of unique kpoints. Append length of kpoint list + 1 for loop below
    kpoint_zero = vcat(findall(!iszero, kpoint_repeating),length(kpoint_repeating)+1)
    for k in 1:length(kpoint_zero)-1
        # Finds overlap of plane waves at unique kpoints
        range = vcat(kpoint_zero[k]:1:(kpoint_zero[k+1])-1)
        S[range,range] = occ_coeff[:,range]'*occ_coeff[:,range] 
    end
    S = (S+S')*0.5
    return S
end


"""
    generate_H(super::Supercell, ehtparams::ehtParams) -> H::Matrix{Float64}

Fills diagonal of the Atomic Orbital Hamiltonian with energy from eht_params.
"""
function generate_H(super::Supercell, ehtparams::ehtParams)
    H = zeros(sum(super.orbitals),sum(super.orbitals))
    for i in 1:length(super.atomlist)
        prev_orb = sum(super.orbitals[1:i-1])
        # d orbitals
        if super.orbitals[i] > 4
            orbital_params = get_eht_params(super.atomlist[i].num, ehtparams)[3]
            for j in 5:9
                H[prev_orb+j,prev_orb+j] = orbital_params.IP
            end
        end
        # p orbitals
        if super.orbitals[i] > 1
            orbital_params = get_eht_params(super.atomlist[i].num, ehtparams)[2]
            for j in 2:4
                H[prev_orb+j,prev_orb+j] = orbital_params.IP
            end
        end
        # s orbital
        orbital_params = get_eht_params(super.atomlist[i].num, ehtparams)[1]
        H[prev_orb+1,prev_orb+1] = orbital_params.IP
    end
    return H
end

"""
    reconstruct_targets_DFT()

The guts of DFTraMO.
"""
# For now, will not use the raMOSystemStatus struct. This will be a to-do to clean up the code!
function reconstruct_targets_DFT(
    psi_target::Array{<:Real},
    num_electrons_left::Int,
    run_name::AbstractString,
    super::Supercell,
    ehtparams::ehtParams,
    wavefxn::ReciprocalWavefunction{3,Float32},
    occ_coeff::Matrix{ComplexF64},
    kpoint_repeating::Vector{Bool},
    G::Vector{Vector{Int64}},
    kptlist::Vector{Vector{Float64}},
    psi_previous::Array{Float64, 3},
    S_original::Matrix{ComplexF64},
    use_prev::Bool,
    prev_mat::AbstractString="",
    )
    
    # Single target run or multiple targets
    num_targets = 1
    if length(size(psi_target)) == 2
        num_targets = size(psi_target)[2]
    end
    num_electrons_left = num_electrons_left - 2*num_targets
    
    #output_name = @sprintf("%s_%05i.mat",run_name,num_electrons_left)
    if use_prev
        # TO DO: write code to load in previous matrix
        open(prev_mat)
        println("Found and loaded remainder file ", prev_mat, ".")
    end
    
    # Calculate overlap between each atomic orbital and planewave
    H = generate_H(super,ehtparams)
    # Checks for spin states
    num_spin_states = size(wavefxn.waves)[1]
    if num_spin_states == 1
        num_spin_up = size(occ_coeff)[2]
        num_spin_down = 0;
    else
        # To do: fix this part once spin up/down is dealt with
        spin_up_coeff = occ_coeff[:,1:num_spin_up]
        spin_down_coeff = occ_coeff[:,1+num_spin_up:num_spin_up+num_spin_down]
    end

    current_orb = 1;
    overlap_target_occupied = zeros(ComplexF64, max(num_spin_up,num_spin_down), num_targets, num_spin_states)

    # Not sure what's going on here just yet; still parsing matlab code (~ ln 852)
    E_mat = Vector{Float64}(undef,0)
    for i in 1:num_targets
        for j in 1:length(super.atomlist)
            prev_orb = sum(super.orbitals[1:j-1])
            # If psi_target is not empty for this atom
            if norm(psi_target[prev_orb+1:prev_orb+super.orbitals[j],i]) > 0
                # Calculate overlap
                overlap_target_temp = calculate_overlap(
                    num_spin_states,
                    num_spin_up,
                    num_spin_down,
                    super.orbitals[j], #num_target_orbitals,
                    occ_coeff,
                    kpoint_repeating,
                    super.atomlist[j].pos, #atom_pos_fract::Vector{Float64},
                    wavefxn.rlatt, #reciprocal_lattice::ReciprocalBasis{3},
                    G,
                    kptlist,#::Vector{Vector{Float64}},
                    get_eht_params(super.atomlist[j].num, ehtparams)
                    )
                    for n in 1:super.orbitals[j]
                        overlap_target_occupied[:,current_orb,:] = overlap_target_occupied[:,current_orb,:]+psi_target[prev_orb+n,i]*overlap_target_temp[:,n,:]
                    end
            end
        end
        E_mat = append!(E_mat, psi_target[:,i]'*H*psi_target[:,i])
        current_orb = current_orb+1
    end
    # Generates correct H and S matrices based on the nature of the analysis (number of spin states, number of remainder states)
    total_num_target_orbs = current_orb-1
    targets_reconstructed = size(psi_previous)[1]-size(psi_previous)[2] # targets already reconstructed
    if num_spin_states == 2
        # spin up
        H = psi_previous[1:num_spin_up,1:(num_spin_up-targets_reconstructed),1]'*overlap_target_occupied[1:num_spin_up,:,1]*diagm(E_mat)*overlap_target_occupied[1:num_spinup,:,1]'*psi_previous[1:num_spin_up,1:(num_spin_up-targets_reconstructed),1]
        S = psi_previous[1:num_spin_up,1:(num_spin_up-targets_reconstructed),1]'*S_original[1:num_spin_up,1:num_spin_up]*psi_previous[1:num_spin_up,1:(num_spin_up-targets_reconstructed),1]
        H = 0.5*(H+H')
        S = 0.5*(S+S')
        (psivect, e_up) = eigen(H,S)
        # sets up new psi_previous or the final raMOs
        tempup = psi_previous[1:num_spin_up,1:(num_spin_up-targets_reconstructed)]*psivect
        # spin down
        H = psi_previous[1:num_spin_down,1:(num_spin_down-targets_reconstructed),2]'*overlap_target_occupied[1:num_spin_down,:,2]*diagm(E_mat)*overlap_target_occupied[1:num_spinup,:,2]'*psi_previous[1:num_spin_down,1:(num_spin_down-targets_reconstructed),2]
        S = psi_previous[1:num_spin_down,1:(num_spin_down-targets_reconstructed),2]'*S_original[1:num_spin_down,1:num_spin_down]*psi_previous[1:num_spin_down,1:(num_spin_down-targets_reconstructed),2]
        H = 0.5*(H+H')
        S = 0.5*(S+S')
        (psivect, e_down) = eigen(H,S)
        tempdown = psi_previous[1:num_spin_down,1:(num_spin_down-targets_reconstructed)]*psivect
        if num_spin_down - targets_reconstructed < 1
            tempdown = tempup+0.0
            e_down = e_up*0.0
        end
    else
        H = psi_previous[:,:,1]'*overlap_target_occupied[:,:,1]*diagm(E_mat)*overlap_target_occupied[:,:,1]'*psi_previous[:,:,1]
        S = psi_previous[:,:,1]'*S_original*psi_previous[:,:,1]
        H = 0.5*(H+H')
        S = 0.5*(S+S')
        (e_up, psi_vect) = eigen(H,S)
        tempup = psi_previous[:,:,1]*psi_vect
    end
    new_remainders = size(psi_previous)[2]-num_targets
    new_psi_previous = zeros(ComplexF64, max(num_spin_up, num_spin_down),new_remainders,num_spin_states);
    if num_spin_up-targets_reconstructed-num_targets > 0
        new_psi_previous[1:num_spin_up,1:(num_spin_up-targets_reconstructed-num_targets),1] = tempup[:,num_targets+1:end]
    end
    print("Number of spin up states remaining: ", num_spin_up-targets_reconstructed-num_targets, "\n")
    psi_up = tempup[:,1:total_num_target_orbs]
    if num_spin_states == 2
        if (num_spin_down-targets_reconstructed-num_targets) > 0
            new_psi_previous[1:num_spin_down,1:(num_spin_down-targets_reconstructed-num_targets),2] = tempup[:,num_targets+1:end]
        end
        print("Number of spin down states remaining: ", num_spin_down-targets_reconstructed-num_targets, "\n")
        psi_down = tempdown[:,1:total_num_target_orbs]
    end
    #write_to_XSF()
    f = open(string(run_name,".xsf"),"w")
    writeXSF(f, data)
    return (new_psi_previous, num_electrons_left)
end

function calculate_overlap(
    num_spin_states::Int,
    num_spin_up::Int,
    num_spin_down::Int,
    num_target_orbitals::Int,
    occ_coeff,
    kpoint_repeating::Vector{Bool},
    atom_pos_fract::StaticArraysCore.SVector{3, Float64},
    reciprocal_lattice::ReciprocalBasis{3},
    G::Vector{Vector{Int64}},
    kptlist::Vector{Vector{Float64}},
    e::Vector{OrbitalParams},
    )
    (num_planewaves, num_occ_states) = size(occ_coeff)
    # Initialize output matrix
    overlap_target_occupied = zeros(num_spin_states, max(num_spin_up, num_spin_down), num_target_orbitals)
    # Initialize secondary overlap container
    overlap = zeros(ComplexF64, num_planewaves*num_occ_states, num_target_orbitals)
    for i in 1:num_occ_states
        if kpoint_repeating[i] == 0 # If this has been done before, copy values.
            overlap[(i-1)*num_planewaves+1:i*num_planewaves,:] = overlap[(i-2)*num_planewaves+1:(i-1)*num_planewaves,:]
        else
            for j in 1:num_planewaves
                direction = -(G[j]+kptlist[i])
                # scaling factor is necessary, bc the STO is not necessarily at the origin
                scalingfactor = exp(2*im*pi*dot(direction, atom_pos_fract))
                direction2 = direction'*reciprocal_lattice
                k_G = norm(direction2) # |k+G|
                if num_target_orbitals >= 1
                    # Deal with s orbitals
                    if k_G == 0 
                        overlap[(i-1)*num_planewaves+j,1] = 0
                    else
                        z = e[1].exp1 # s orbital, First zeta value
                        # HARDCODED INTEGRALS OF radial_STOs x spherical bessel function j0(k_G*r)
                        if e[1].n_quant == 1
                            overlap[(i-1)*num_planewaves+j,1] = scalingfactor * N_L(0) * 4*z^(5/2)/(k_G^2 + z^2)^2
                        elseif e[1].n_quant == 2
                            overlap[(i-1)*num_planewaves+j,1] = scalingfactor * N_L(0) * -(4*z^(5/2)*(k_G^2 - 3*z^2))/(sqrt(3) * (k_G^2 + z^2)^3)
                        elseif e[1].n_quant == 3
                            overlap[(i-1)*num_planewaves+j,1] = scalingfactor * N_L(0) * (16*sqrt(2/5)*z^(9/2)*(z^2 - k_G^2))/(k_G^2 + z^2)^4
                        elseif e[1].n_quant == 4
                            overlap[(i-1)*num_planewaves+j,1] = scalingfactor * N_L(0) * (16*z^(9/2)*(k_G^4 - 10*k_G^2*z^2 + 5*z^4))/(sqrt(35)*(k_G^2 + z^2)^5)
                        elseif e[1].n_quant == 5
                            overlap[(i-1)*num_planewaves+j,1] = scalingfactor * N_L(0) * (32*sqrt(2/7)*z^(13/2)*(3*k_G^4 - 10*k_G^2*z^2 + 3*z^4))/(3*(k_G^2 + z^2)^6)
                        elseif e[1].n_quant == 6
                            overlap[(i-1)*num_planewaves+j,1] = scalingfactor * N_L(0) * -(32*sqrt(2/231)*z^(13/2)*(k_G^6 - 21*k_G^4*z^2 + 35*k_G^2*z^4 - 7*z^6))/(k_G^2 + z^2)^7
                        elseif e[1].n_quant == 7
                            overlap[(i-1)*num_planewaves+j,1] = scalingfactor * N_L(0) * (512*z^(17/2)*(-k_G^6 + 7*k_G^4*z^2 - 7*k_G^2*z^4 + z^6))/(sqrt(429)*(k_G^2 + z^2)^8)
                        end
                    end
                end
                if num_target_orbitals >= 4
                    # Deal with p orbitals. Store overlap in last p orbital for now.
                    if k_G == 0 
                        overlap[(i-1)*num_planewaves+j,2:4] .= 0
                    else
                        z = e[2].exp1 # p orbital, First zeta value
                        if e[2].n_quant == 2
                            overlap[(i-1)*num_planewaves+j,4] = scalingfactor * N_L(1) * (16*k_G*z^(7/2))/(sqrt(3)*(k_G^2 + z^2)^3) # for some reason radial_overlap()'s function does not give the correct answer
                        elseif e[2].n_quant == 3
                            overlap[(i-1)*num_planewaves+j,4] = scalingfactor * N_L(1) * -(16*sqrt(2/5)*k_G*z^(7/2)*(k_G^2 - 5*z^2))/(3*(k_G^2 + z^2)^4)
                        elseif e[2].n_quant == 4
                            overlap[(i-1)*num_planewaves+j,4] = scalingfactor * N_L(1) * (32*k_G*z^(11/2)*(-3*k_G^2 + 5*z^2))/(sqrt(35)*(k_G^2 + z^2)^5)
                        elseif e[2].n_quant == 5
                            overlap[(i-1)*num_planewaves+j,4] = scalingfactor * N_L(1) * (32*sqrt(2/7)*k_G*z^(11/2)*(3*k_G^4 - 42*k_G^2*z^2 + 35*z^4))/(15*(k_G^2 + z^2)^6)
                        elseif e[2].n_quant == 6
                            overlap[(i-1)*num_planewaves+j,4] = scalingfactor * N_L(1) * (256*sqrt(2/231)*k_G*z^(15/2)*(3*k_G^4 - 14*k_G^2*z^2 + 7*z^4))/(3*(k_G^2 + z^2)^7)
                        end
                        # Overlap is calculated for pz orbital along direction of planewave.
                        # Now determine overlap of px, py, pz orbitals in reference frame
                        # Φ = angle between projection of k+G on xy-plane to +x-axis
                        # Θ = angle between k+G to +z-axis
                        if sqrt(direction2[1]^2+direction2[2]^2) < 0.00001 # If vector lies along z direction
                            cosΦ = 1; sinΦ = 0; sinΘ = 0;
                        else
                            sinΦ = direction2[2]/sqrt(direction2[1]^2+direction2[2]^2)
                            cosΦ = direction2[1]/sqrt(direction2[1]^2+direction2[2]^2)
                            sinΘ = sqrt(direction2[1]^2+direction[2]^2)/k_G
                        end
                        cosΘ = direction2[3]/k_G
                        # Rotation matrix for p orbitals
                        C_1 = [sinΘ*cosΦ, sinΘ*sinΦ, cosΘ]
                        # Apply rotation to p orbitals
                        overlap[(i-1)*num_planewaves+j,2:4] = overlap[(i-1)*num_planewaves+j,4].*C_1
                    end
                end
                if num_target_orbitals >= 9
                    # Deal with d orbitals. Store overlap in last d orbital for now.
                    if k_G == 0 
                        overlap[(i-1)*num_planewaves+j,5:9] .= 0
                    else
                        z1 = e[3].exp1 # d orbital, First zeta value.
                        z2 = e[3].exp2 # d orbital, Second zeta value.
                        c1 = e[3].coeff1
                        c2 = e[3].coeff2
                        dznorm = ((4*z1*z2)/(z1 + z2)^2)^(e[3].n_quant+1/2) # double-zeta normalization factor
                        # scaling_factor * N_L * double-zeta_normalization_factor * integral_btwn_double_zeta_STO_and_spherical_bessel_j3()
                        if e[3].n_quant == 3
                            overlap[(i-1)*num_planewaves+j,9] = scalingfactor * N_L(2) * dznorm * (32*sqrt(2/5)*k_G^2*(c2*(k_G^2 + z1^2)^4*z2^(9/2) + c1*z1^(9/2)*(k_G^2 + z2^2)^4))/((k_G^2 + z1^2)^4 * (k_G^2 + z2^2)^4)
                        elseif e[3].n_quant == 4
                            overlap[(i-1)*num_planewaves+j,9] = scalingfactor * N_L(2) * dznorm * -(32/sqrt(35)*k_G^2*((c1*z1^(9/2)*(k_G^2 - 7*z1^2))/(k_G^2 + z1^2)^5 + (c2*z2^(9/2)*(k_G^2 - 7*z2^2))/(k_G^2 + z2^2)^5))
                        elseif e[3].n_quant == 5
                            overlap[(i-1)*num_planewaves+j,9] = scalingfactor * N_L(2) * dznorm * -(256/15*sqrt(2/7)*k_G^2*((c1*z1^(13/2)*(3*k_G^2 - 7*z1^2))/(k_G^2 + z1^2)^6 + (c2*z2^(13/2)*(3*k_G^2 - 7*z2^2))/(k_G^2 + z2^2)^6))
                        end
                        if sqrt(direction2[1]^2+direction2[2]^2)<0.00001 # if the vector lies along the z direction
                            cosΦ = 1; sinΦ = 0; sinΘ = 0;
                        else
                            sinΦ = direction2[2]/sqrt(direction2[1]^2+direction[2]^2)
                            cosΦ = direction2[1]/sqrt(direction2[1]^2+direction[2]^2)
                            sinΘ = sqrt(direction2[1]^2+direction[2]^2)/k_G
                        end
                        cosΘ = direction2[3]/k_G
                        # d orbital rotation matrix
                        C_2 = [
                            sqrt(3)/2*sinΘ^2*(cosΦ^2 - sinΦ^2),
                            1 - 3/2*sinΘ*2,
                            sqrt(3)*cosΦ*sinΘ^3,
                            sqrt(3)*cosΘ*sinΘ*cosΦ,
                            sqrt(3)*cosΘ*sinΘ*sinΦ
                            ]
                        # apply rotation matrix
                        overlap[(i-1)*num_planewaves+j,5:9] = overlap[(i-1)*num_planewaves+j,9].*C_2
                    end
                end
            end
        end
    end
    # Create overlap matrix
    overlap_target_occupied = zeros(ComplexF64, num_occ_states, num_target_orbitals)
    if num_spin_states == 1
        for i in 1:num_occ_states
            overlap_target_occupied[i,1:num_target_orbitals] = occ_coeff[:,i]'*overlap[(i-1)*num_planewaves+1:i*num_planewaves,:]
        end
    else # spin states = 2
        for i in 1:num_spin_up
            overlap_target_occupied[i,:,1] = spin_up_coeff[:,i]'overlap[(i-1)*num_planewaves+1:i*num_planewaves,:]
        end
        for i in num_spin_up+1:num_spin_up+num_spin_down
            overlap_target_occupied[i-num_spin_up,:,2] = spin_down_coeff[:,i-num_spin_up]'*overlap[(i-1)*num_planewaves+1:i*num_planewaves,:]
        end
    end
return overlap_target_occupied
end
"""
    N_L(l::Int)]

Returns a constant of the planewave expansion, equal to (im^l)*(4*pi*(2l+1))^0.5.
"""
function N_L(l::Int)
    return (im^l)*(4*pi*(2l+1))^0.5
end