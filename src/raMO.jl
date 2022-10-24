"""
get_eht_params(atom_num::Integer, atom_orbital::Integer, eht_params::ehtParams) -> OrbitalParams

Search for parameters in the loaded ehtParams for corresponding atom. atom_orbital corresponds to
the quantum number l.
"""
function get_eht_params(atom_num::Integer, atom_orbital::Integer, eht_params::ehtParams)
    return eht_params.data[atom_num, atom_orbital+1]
end

"""
    make_overlap_mat(occ_coeff, kpoint_repeating::AbstractVector{Bool})

Creates overlap matrix with occupied coefficients. If the kpoints are the same, they can overlap.
Otherwise, they do not.
"""
function make_overlap_mat(occ_coeff, kpoint_repeating::AbstractVector{Bool})
    num_occ_states = length(kpoint_repeating)
    S = zeros(num_occ_states,num_occ_states)
    kpoint_zero = findall(x -> x == 0, kpoint_repeating)
    S[kpoint_zero,kpoint_zero] = occ_coeff[:,kpoint_zero]'*occ_coeff[:,kpoint_zero]
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
            orbital_params = get_eht_params(super.atomlist[i].num, 2, ehtparams)
            for j in 5:9
                H[prev_orb+j,prev_orb+j] = orbital_params.IP
            end
        end
        # p orbitals
        if super.orbitals[i] > 1
            orbital_params = get_eht_params(super.atomlist[i].num, 1, ehtparams)
            for j in 2:4
                H[prev_orb+j,prev_orb+j] = orbital_params.IP
            end
        end
        # s orbital
        orbital_params = get_eht_params(super.atomlist[i].num, 0, ehtparams)
        H[prev_orb+1,prev_orb+1] = orbital_params.IP
    end
    return H
end

"""
    reconstruct_targets_DFT(
        psi_target::AbstractArray{<:Real},
        num_electrons_left::Integer,
        run_name::AbstractString,
        super::Supercell,
        ehtparams::ehtParams,
        wavefxn::ReciprocalWavefunction{3,Float32},
        use_prev::Bool,
        prev_mat::AbstractString=""
    )

The guts of DFTraMO.
"""
# For now, will not use the raMOSystemStatus struct. This will be a to-do to clean up the code!
function reconstruct_targets_DFT(
    psi_target::AbstractArray{<:Real},
    num_electrons_left::Integer,
    run_name::AbstractString,
    super::Supercell,
    ehtparams::ehtParams,
    wavefxn::ReciprocalWavefunction{3,Float32},
    use_prev::Bool,
    prev_mat::AbstractString=""
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
        open(prev_mat,"r")
        println("Found and loaded remainder file ", prev_mat, ".")
    end
    
    # Calculate overlap between each atomic orbital and planewave
    H = generate_H(super,ehtparams)
    # Checks for spin states
    num_spin_states = size(wavefxn.waves)[1]
    if num_spin_states == 1
        spin_up_coeff = 1
        spin_down_coeff = 1
    else
        # To do: fix this part once spin up/down is dealt with
        spin_up_coeff = occ_coeff[:,1:num_spin_up]
        spin_down_coeff = occ_coeff[:,1+num_spin_up:num_spin_up+num_spin_down]
    end

    # Not sure what's going on here just yet; still parsing matlab code (~ ln 852)
    for i in 1:num_targets
        for j in 1:length(super.atomlist)
            prev_orb = sum(super.orbitals[1:j-1])
            # If psi_target is not empty for this atom
            if norm(psi_target[prev_orb+1:prev_orb+super.orbitals[j],i]) > 0
                @info norm(psi_target[prev_orb+1:prev_orb+super.orbitals[j],i])
                # Calculate overlap
                #== overlap_target_temp, E_mat_temp = calculate_overlap(
                    num_spin_states,
                    num_spin_up,
                    num_spin down,
                    super.orbitals[j]), # number of target orbitals for the atom
                    num_occ_states,
                    kpoint_repeating,
                    super.atomlist[j].pos, # atomic position in fractional coordinates
                    super.waves[1].basis # reciprocal lattice vectors
                    )==#
            end
        end
    end
end

# Add a docstring here...
function calculate_overlap(
    num_spin_states::Integer,
    num_spin_up::Integer,
    num_spin_down::Integer,
    num_target_orbitals::Integer,
    num_occ_states::Integer,
    kpoint_repeating::AbstractVector{Bool},
    atom_pos_fract::AbstractVector{<:Real},
    reciprocal_lattice::ReciprocalBasis{3},
    G::AbstractVector{<:Integer},
    kptlist::AbstractVector{<:AbstractVector{<:Real}},
    e::AbstractVector{OrbitalParams},
)
    # Initialize output matrix
    overlap_target_occupied = zeros(num_spin_states, max(num_spin_up, num_spin_down), num_target_orbitals)
    # Initialize secondary overlap container
    overlap = zeros(num_planewaves*num_occ_states, num_target_orbitals)
    for i in 1:num_occupied_states
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
                    # Deal with p orbitals. Store overlap in pz orbital for now.
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
                    end
                end
                if num_target_orbitals >= 9
                    # Deal with d orbitals. Store overlap in dz2 orbital for now.
                    if k_G == 0 
                        overlap[(i-1)*num_planewaves+j,5:9] .= 0
                    end
                    z1 = e[3].exp1 # d orbital, First zeta value.
                    z2 = e[3].exp2 # d orbital, Second zeta value.
                    c1 = e[3].coeff1
                    c2 = e[3].coeff2
                    if e[3].n_quant == 3
                        overlap[(i-1)*num_planewaves+j,6] = scalingfactor * N_L(2) * (32*sqrt(2/5)*k_G^2*(c2*(k_G^2 + z1^2)^4*z2^(9/2) + c1*z1^(9/2)*(k_G^2 + z2^2)^4))/((k_G^2 + z1^2)^4 * (k_G^2 + z2^2)^4)
                    end
                end
            end
        end
    end
end
"""
    N_L(l::Integer)

Returns a constant of the planewave expansion, equal to `(im^l)*(4*pi*(2l+1))^0.5`.
"""
function N_L(l::Integer)
    return (im^l)*(4*pi*(2l+1))^0.5
end