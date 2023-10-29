# TODO: is there a reason why neither orb nor val_e are constants/have no docstring?
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

const ORB_DICT = Dict{String, Int}(
    [Electrum.ELEMENTS[n] => orb[n] for n in eachindex(orb)]
)

const E_DICT = Dict{String, Int}(
    [Electrum.ELEMENTS[n] => val_e[n] for n in eachindex(val_e)]
)

# TODO: can we use regex to handle this more efficiently?
const AO_RUNS = Dict(
    "s"     => 1,
    "px"    => 2,
    "py"    => 3,
    "pz"    => 4,
    "dx2y2" => 5,
    "dx2-y2"=> 5,
    "dz2"   => 6,
    "dxy"   => 7,
    "dxz"   => 8,
    "dyz"   => 9
)

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
    return new(atom_num, valence, l_quant, n_quant, IP,exp1, exp2, coeff1, coeff2)
    end
end

Base.getindex(o::OrbitalParams, s::Symbol) = getproperty(o,s)

"""
    ehtParams <: AbstractMatrix{OrbitalParams}

Matrix of OrbitalParams (atom x ang).
"""
struct ehtParams <: AbstractMatrix{OrbitalParams}
    data::Matrix{OrbitalParams}
    ehtParams(data::AbstractMatrix{OrbitalParams}) = new(data)
end

Base.size(e::ehtParams) = size(e.data)
Base.axes(e::ehtParams) = axes(e.data)
Base.getindex(e::ehtParams, i...) = e.data[i...]

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

Electrum.basis(s::Supercell) = basis(s.atomlist)
Electrum.PeriodicAtomList(s::Supercell) = s.atomlist
Base.getindex(s::Supercell, i...) = getindex(s.atomlist, i...)

"""
    OccupiedStates(
        coeff::AbstractMatrix{<:Number},
        kpt::AbstractMatrix{<:AbstractVector},
        G::AbstractMatrix{<:AbstractVector}
    )

Returns an `OccupiedStates` struct. The coeff matrix is `num_occupied_states` by `num_occupied_pw`
in dimensions, while kpt is `num_occupied_states` in length, and G is `num_occupied_pw` in length.
"""
struct OccupiedStates
    coeff::Matrix{ComplexF32}
    kpt::Matrix{SVector{3, Float64}}
    G::Matrix{SVector{3, Int64}}
    function OccupiedStates(
        coeff::AbstractMatrix{<:Number},
        kpt::AbstractMatrix{<:AbstractVector},
        G::AbstractMatrix{<:AbstractVector}
    )
        return new(coeff, kpt, G)
    end
end

"""
    DFTraMO.RunInfo

Contains information from a run entry in the YAML input file. All runs in the file are collected
a single `Vector{RunInfo}` object when `dftramo_run()` is called.
"""
struct RunInfo
    name::String
    type::String
    site_file::String
    sites::Vector{Int}
    radius::Float64
    rsphere::Float64
    function RunInfo(
        name::AbstractString,
        type::AbstractString,
        site_file::AbstractString,
        sites::AbstractVector{<:Integer},
        radius::Number,
        rsphere::Number
    )
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