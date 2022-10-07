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

function read_eht_params()
    ln = readlines(open("DFT_raMO_eht_parms.dat","r"))
    ln = ln[4:length(ln)]
    ln = filter(!isempty,split.(ln))
    unique_atoms_index = unique(i->mapreduce(permutedims,vcat,ln)[:,1][i],eachindex(mapreduce(permutedims,vcat,ln)[:,1]))
    mat = ehtParams(Matrix{OrbitalParams}(undef,length(ln),4))
    atom_counter = 0
    for i in 1:length(ln)
        current = ln[i]
        # Assigns atom_num,valence,l_quant,n_quant,IP,exp1,exp2,coeff1,coeff2
        orbs = Dict("s"=>0,"p"=>1,"d"=>2,"f"=>3)
        orb_param = OrbitalParams(parse(Int,current[2]), parse(Int,current[3]), get(orbs,current[6],0), parse(Int,current[5]), parse(Float64,current[7]), parse(Float64,current[8]), parse(Float64,current[9]), parse(Float64,current[10]), parse(Float64,current[11]))
        if i in unique_atoms_index
            atom_counter = atom_counter + 1
        end
        mat.data[atom_counter, get(orbs,current[6],0)+1] = orb_param
    end
    return mat
end