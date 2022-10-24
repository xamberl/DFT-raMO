"""
    OrbitalParams

Stores parameters for Slater-type orbitals.
"""
struct OrbitalParams
    atom_num::Int
    valence::Int
    l_quant::Int
    n_quant::Int
    IP::Float64
    exp1::Float64
    exp2::Float64
    coeff1::Float64
    coeff2::Float64
    function OrbitalParams(
        atom_num::Integer,
        valence::Integer,
        l_quant::Integer,
        n_quant::Integer,
        IP::Real,
        exp1::Real,
        exp2::Real,
        coeff1::Real,
        coeff2::Real
    )
    return new(atom_num,valence,l_quant,n_quant,IP,exp1,exp2,coeff1,coeff2)
    end
end

function Base.getindex(o::OrbitalParams, s::Symbol)
    return getproperty(o,s)
end

"""
    ehtParams

Matrix of OrbitalParams (atom x ang)
"""
struct ehtParams
    data::Matrix{OrbitalParams}
    function ehtParams(data::Matrix{OrbitalParams})
        return new(data)
    end
end

function Base.getindex(e::ehtParams, atom::Integer, ang::Integer)
    return e.data[atom, ang]
end

"""
    supercell(atomlist::Atomlist{3}, orblist_by_type::Dict{String,Int64})

Returns a supercell struct AtomList{3} and a Vector{Int} with number of orbitals associated with each atom
"""
struct Supercell
    atomlist::AtomList{3}
    orbitals::Vector{Int}
    function Supercell(atomlist::AtomList{3}, orblist_by_type::Dict{String,Int64})
        orbitals = Vector{Int}(undef,0)
        for atom in atomlist.coord
            if !haskey(orblist_by_type, atom.name)
                error("Number of orbitals for atom", atom.name,"cannot be found.")
            end
            # Get the number of orbitals for current atom and add to sum
            push!(orbitals,get(orblist_by_type,atom.name,0))
        end
        return new(atomlist,orbitals)
    end
end

#=
Might want this to make reconstruct_targets_DFT() neater?
"""
"""
struct raMOSystemStatus
    num_electrons_left::Integer
    num_occ_states::Integer
    num_planewaves::Integer
    kptlist::KPointList{3}
    G::AbstractVector{<:Integer}
    #occ_coeff::
end
=#
