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
        push!(at_lab,i[1])
    end
    test = findall(x->x=="nothing",at_lab)
    if isempty(test)
        error("Could not find element ", atom_symbol, " in DFT_raMO_eht_parms.dat file.")
    else
    end
end==#

function read_eht_params()
    ln = readlines(open("DFT_raMO_eht_parms.dat","r"))
    ln = ln[4:length(ln)]
    ln = filter(!isempty,split.(ln))
end