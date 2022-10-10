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