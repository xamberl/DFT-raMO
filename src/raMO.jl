"""
    get_eht_params(atom_num::Int) -> OrbitalParams

Search for parameters in the loaded ehtParams for corresponding atom.
"""
function get_eht_params(atom_num::Int, atom_orbital::Int, eht_params::ehtParams)
    return eht_params.data[atom_num, atom_orbital]
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