"""
get_eht_params(atom_num::Int, atom_orbital::Int, eht_params::ehtParams) -> OrbitalParams

Search for parameters in the loaded ehtParams for corresponding atom. atom_orbital corresponds to the quantum number l.
"""
function get_eht_params(atom_num::Int, atom_orbital::Int, eht_params::ehtParams)
    return eht_params.data[atom_num, atom_orbital+1]
end


"""
    make_overlap_mat(occ_coeff, kpoint_repeating::Vector{Bool})

    Creates overlap matrix with occupied coefficients
    If the kpoints are the same, they can overlap. Otherwise, they do not.
"""
function make_overlap_mat(occ_coeff, kpoint_repeating::Vector{Bool})
    num_occ_states = length(kpoint_repeating)
    S = zeros(num_occ_states,num_occ_states)
    kpoint_zero = findall(x -> x == 0, kpoint_repeating)
    S[kpoint_zero,kpoint_zero] = occ_coeff[:,kpoint_zero]'*occ_coeff[:,kpoint_zero]
    return S
end

"""
    make_target_AO(make_target_AO(atom_site::Vector{Int}, orbital_to_use::Int, total_num_orbitals::Int) -> Vector

Returns a vector of length total_num_orbitals with "1" in the corresponding atomic orbital
"""
#==
In the original DFTraMO, there is a make_target_massAO function, where instead of returning a
vector with only one "1" in psi_target, it has multiple "1"s in psi_target. This is resolved by
making atom_site a Vector, so we can get it for one site (length(atom_site) ==  1) or multiple.
==#
function make_target_AO(atom_site::Vector{Int}, target_orbital::Int, super::Supercell)
    # psi_target has a length of total number of orbitals in the supercell
    psi_target = zeros(sum(super.orbitals))
    for atom in atom_site
        psi_target[sum(super.orbitals[1:atom])-super.orbitals[atom]+target_orbital] = 1
    end
    return psi_target
end

"""
    make_target_cluster_sp(site_list::Vector{Vector{Real}}, radius::Real, site_num::Int, super::Supercell)

Returns a vector of length total_num_orbitals with "1" in the corresponding s & p atomic orbitals
if they are within specified radius to the void.
"""
function make_target_cluster_sp(site_list::Vector{Vector{Float64}}, radius::Real, site_num::Int, super::Supercell)
    # psi_target has a length of total number of orbitals in the supercell
    psi_target = zeros(sum(super.orbitals))
    # Loops through every atom and checks to see if it's within the radius to the void site
    for n in 1:length(super.atomlist)
        for j in -1:2
            for k in -1:2
                for l in -1:2
                    # Translation of atom in cartesian coordinates
                    new_pos = super.atomlist.basis*super.atomlist[n].pos .+ super.atomlist.basis*[j,k,l]
                    check_distance = norm(new_pos-site_list[site_num])
                    if check_distance <= radius
                        # Δr is used to weigh the p orbitals in each direction
                        # Is there a reason why this is negative?
                        Δr = -(new_pos-site_list[site_num])/check_distance
                        # Fill in corresponding s orbital
                        psi_target[sum(super.orbitals[1:n])-super.orbitals[n]+1] = 0.5^0.5
                        # Fill in corresponding p orbitals, if any
                        if (super.orbitals[n] > 1)
                            for i in 1:3
                                psi_target[sum(super.orbitals[1:n])-super.orbitals[n]+1+i] = 0.5^(0.5*Δr[i])
                            end
                        end
                    end
                end
            end
        end
    end
    return psi_target
end

"""
    make_target_hybrid()

Makes a custom hybrid with manually specified coefficients from a cluster file that corresponds to an atom
"""
#==
Currently echoes MATLAB version of DFT-raMO. Cluster file should have one Float per line,
in multiples of 9.
==#
function make_target_hybrid(cluster_list::Vector{Vector{Float64}}, atom_num::Int, super::Supercell)
    num_targets = length(cluster_list)/9
    psi_target = zeros(sum(super.orbitals), num_targets)
    AO_number = sum(super.orbitals[1:atom_num])-super.orbitals[atom_num]
    # Loops through targets
    for n in 1:num_targets
        # Assigns s coefficient
        psi_target[AO_number+1,n] = cluster_list[(n-1)*9+1]
        # Assigns p coefficients
        if super.orbitals[atom_num] > 1
            for p in 1:3
                psi_target[AO_number+1+p,n] = cluster_list[(n-1)*9+1+p]
            end
        end
        # Assigns d coefficients
        if super.orbitals[atom_num] > 4
            for d in 1:3
                psi_target[AO_number+4+d,n] = cluster_list[(n-1)*9+4+d]
            end
        end
    end
    return psi_target
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

function calculate_overlap(
    num_spin_states::Int,
    num_spin_up::Int,
    num_spin_down::Int,
    num_target_orbitals::Int,
    num_occ_states::Int,
    kpoint_repeating::Vector{Bool},
    atom_pos_fract::Vector{Float64},
    reciprocal_lattice::ReciprocalBasis{3},
    G::Vector{Int},
    kptlist::Vector{Vector{Float64}},
    e::OrbitalParams,
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
                # Vince's method
                direction = -(G[j]+kptlist[i])
                scalingfactor = exp(2*im*pi*dot(direction, atom_pos_fract))
                direction2 = direction'*reciprocal_lattice
                direction_norm = norm(direction2)
                if num_target_orbitals > 1
                    # Deal with s orbitals
                    if direction_norm == 0 
                        overlap[(i-1)*num_planewaves+j,1] = 0
                    else
                        if e.n_quant == 1
                            overlap[(i-1)*num_planewaves+j,1] = scalingfactor*8*e.exp1^2/(direction_norm^2+e.exp1^2)^2*(e.exp1*pi)^0.5
                        elseif e.n_quant == 2
                        elseif e.n_quant == 3
                        elseif e.n_quant == 4
                            overlap[(i-1)*num_planewaves+j,1] = scalingfactor*96*e.exp1^4*(direction_norm^4-10*direction_norm^2*e.exp1^2+5*e.exp1^4)/(direction_norm^2+e.exp1^2)^5*(pi*e.exp1/315)^0.5
                        elseif e.n_quant == 5
                        elseif e.n_quant == 6
                        elseif e.n_quant == 7
                        end
                    end
                end
            end
        end
    end
end

function N_L(l::Int)
    return (im^l)*(4*pi*(2l+1))^0.5
end