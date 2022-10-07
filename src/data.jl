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
        atom_num::Int,
        valence::Int,
        l_quant::Int,
        n_quant::Int,
        IP::Float64,
        exp1::Float64,
        exp2::Float64,
        coeff1::Float64,
        coeff2::Float64
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

