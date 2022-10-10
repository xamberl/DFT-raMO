"""
    get_eht_params(atom_num::Int) -> OrbitalParams

Search for parameters in the loaded ehtParams for corresponding atom.
"""
function get_eht_params(atom_num::Int, atom_orbital::Int, eht_params::ehtParams)
    return eht_params.data[atom_num, atom_orbital]
end

# Creates overlap matrix with occupied coefficients
# If the kpoints are the same, they can overlap. Otherwise, they do not.
function make_overlap_mat(occ_coeff, kpoint_repeating::Vector{Bool})
    num_occ_states = length(kpoint_repeating)
    S = zeros(num_occ_states,num_occ_states)
    kpoint_zero = findall(x -> x == 0, kpoint_repeating)
    S[kpoint_zero,kpoint_zero] = occ_coeff[:,kpoint_zero]'*occ_coeff[:,kpoint_zero]
    return S
end

"""
    make_target_AO(make_target_AO(atom_site::Int, orbital_to_use::Int, total_num_orbitals::Int) -> Vector

Returns a vector of length total_num_orbitals with "1" in the corresponding atomic orbital
"""
function make_target_AO(atom_site::Int, target_orbital::Int, super::Supercell)
    # psi_target has a length of total number of orbitals in the supercell
    psi_target = zeros(sum(super.orbitals))
    psi_target[sum(super.orbitals[1:atom_site])-super.orbitals[atom_site]+target_orbital] = 1
    return psi_target
end