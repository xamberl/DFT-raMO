#=="""
    get_eht_params

Search for parameters in DFT_raMO_eht_parms.dat file for corresponding atom.
"""
function get_eht_params(atom_symbol::String)
    ln = readlines(open("DFT_raMO_eht_parms.dat","r"))
    ln = ln[4:length(ln)]
    ln = filter(!isempty,split.(ln))
    at_lab = Vector{String}(undef,0)
    for i in ln
        push!(at_lab,current[1])
    end
    test = findall(x->x=="nothing",at_lab)
    if isempty(test)
        error("Could not find element ", atom_symbol, " in DFT_raMO_eht_parms.dat file.")
    else
    end
end==#

# Creates overlap matrix with occupied coefficients
# If the kpoints are the same, they can overlap. Otherwise, they do not.
function make_overlap_mat(occ_coeff, kpoint_repeating::Vector{Bool})
    num_occ_states = length(kpoint_repeating)
    S = zeros(num_occ_states,num_occ_states)
    kpoint_zero = findall(x -> x == 0, kpoint_repeating)
    S[kpoint_zero,kpoint_zero] = occ_coeff[:,kpoint_zero]'*occ_coeff[:,kpoint_zero]
    return S
end