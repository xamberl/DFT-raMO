orb = [
    1, 1, 1, 1, 4, 4, 4, 4, 4, 4,
    1, 1, 4, 4, 4, 4, 4, 4,
    1, 1, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 4, 4, 4, 4, 4, 9,
    1, 1, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 4, 4, 4, 4, 4, 9,
    1, 1, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 4, 4, 4, 4, 4, 4,
    1, 1, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
    16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 
    16, 16, 16, 16, 16, 16, 16, 16
]

val_e = [
    1, 2, 1, 2, 3, 4, 5, 6, 7, 8,
    1, 2, 3, 4, 5, 6, 7, 8,
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 3, 4, 5, 6, 7, 8,
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 3, 4, 5, 6, 7, 8,
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
    4, 5, 6, 7, 8, 9, 10, 11, 12, 3, 4, 5, 6, 7, 8,
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
    4, 5, 6, 7, 8, 9, 10, 11, 12, 3, 4, 5, 6, 7, 8,
]

const orb_dict = Dict{String, Int}(
    [Electrum.ELEMENTS[n] => orb[n] for n in eachindex(orb)]
)

const e_dict = Dict{String, Int}(
    [Electrum.ELEMENTS[n] => val_e[n] for n in eachindex(val_e)]
)

const AO_RUNS = Dict("s"=>1, "px"=>2, "py"=>3, "pz"=>4, "dx2y2"=>5, "dx2-y2"=>5, "dz2"=>6, "dxy"=>7, "dxz"=>8, "dyz"=>9)

const CAGE_RUNS = [
    "sp" # hybrid cage states by distance
]

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

"""
    Supercell(atomlist::Atomlist{3}, orblist_by_type::Dict{String,Int64})

Returns a supercell struct AtomList{3} and a Vector{Int} with number of orbitals associated with each atom
"""
struct Supercell
    atomlist::PeriodicAtomList{3}
    orbitals::Vector{Int}
    function Supercell(atomlist::PeriodicAtomList{3}, orblist_by_type::Dict{String,Int64})
        orbitals = Vector{Int}(undef,0)
        for atom in atomlist.atoms
            if !haskey(orblist_by_type, atom.atom.name)
                error("Number of orbitals for atom", atom.atom.name,"cannot be found.")
            end
            # Get the number of orbitals for current atom and add to sum
            push!(orbitals,get(orblist_by_type,atom.atom.name,0))
        end
        return new(atomlist,orbitals)
    end
end

"""
    OccupiedStates(coeff::Matrix{ComplexF32}, kpt::Vector{SVector{3, Float64}}, G::Vector{Vector{Int64}})

Returns an OccupiedStates struct. The coeff matrix is num_occupied_states by num_occupied_pw in dimensions,
while kpt is num_occupied_states in length, and G is num_occupied_pw in length.
"""
struct OccupiedStates
    coeff::Array{ComplexF32}
    kpt::Array{SVector{3, Float64}}
    G::Array{SVector{3, Int64}}
    function OccupiedStates(coeff::Array{ComplexF32}, kpt::Array{SVector{3, Float64}}, G::Array{SVector{3, Int64}})
        return new(coeff, kpt, G)
    end
end

struct RunInfo
    name::AbstractString
    type::AbstractString
    site_file::AbstractString
    sites::AbstractVector{Int}
    radius::AbstractFloat
    rsphere::AbstractFloat
    function RunInfo(name::AbstractString, type::AbstractString, site_file::AbstractString, sites::AbstractVector{Int}, radius::AbstractFloat, rsphere::AbstractFloat)
        return new(name, type, site_file, sites, radius, rsphere)
    end
end

#==struct raMOStatus
    fermi::NamedTuple{(:fermi, :alphabeta), Tuple{Float64, Float64}}
    xtal::PeriodicAtomList{3}
    geo::PeriodicAtomList{3}
    kpt::AbstractVector{Int}
    wave::Planewavefunction
    num_raMO::Int
    num_electrons_left::Int
end==#