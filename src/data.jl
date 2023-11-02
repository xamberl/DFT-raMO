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
        @assert all(sign.(sites) .> 0) "Some site values are not positive integers."
        @assert sign(radius) >= 0 "Radius is a negative value."
        @assert sign(rsphere) >= 0 "Rsphere is a negative value."
        return new(name, type, site_file, sites, radius, rsphere)
    end
end

Base.size(e::ehtParams) = size(e.data)
Base.axes(e::ehtParams) = axes(e.data)
Base.getindex(e::ehtParams, i...) = e.data[i...]

"""
    raMOCheckpoint{T}

Contains checkpoint data, a 3D array of remainder coefficients, as well as the number of electrons
left and the number of the current raMO in the sequence.
"""
struct raMOCheckpoint{T<:Number} <: AbstractArray{T,3}
    coeff::Array{T,3}
    electrons_left::Int
    num_ramos::Int
    function raMOCheckpoint{T}(
        coeff::AbstractArray{3},
        electrons_left::Integer,
        num_ramos::Integer
    ) where T
        #@assert electrons_left >= 0 "The number of electrons left is a negative value."
        @assert num_ramos >= 0 "The number of raMOs reconstructed is a negative value."
        return new(coeff, electrons_left, num_ramos)
    end
end

function raMOCheckpoint(
    psi::AbstractArray{T,3},
    electrons_left::Integer,
    num_ramos::Integer
) where T
    return raMOCheckpoint{T}(psi, electrons_left, num_ramos)
end

Base.size(x::raMOCheckpoint) = size(x.coeff)
Base.getindex(x::raMOCheckpoint, i...) = x.coeff[i]

"""
    raMODFTData

Contains all of the crystal and wavefunction information needed to perform a DFT-raMO run.
"""
struct raMODFTData
    xtal::Crystal{3}
    wave::PlanewaveWavefunction{3,ComplexF32}
    fermi::Float64
end

Electrum.basis(x::raMODFTData) = basis(x.xtal.atoms)
Electrum.Crystal(x::raMODFTData) = x.xtal
Electrum.PeriodicAtomList(x::raMODFTData) = x.xtal.atoms
Electrum.PlanewaveWavefunction(x::raMODFTData) = x.wave
Electrum.fermi(x::raMODFTData) = x.fermi
kptmesh(x::raMODFTData) = diag(x.xtal.transform)
Electrum.supercell(x::raMODFTData) = supercell(x.xtal.atoms, kptmesh(x))

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
    raMOInput

Contains all of the information supplied by the user in the YAML file with raMO input data. This
includes the DFT data (crystal, wavefunction, Fermi energy), run list, energy ranges for orbital
reconstruction, path to the checkpoint file, and whether to automatically use `Psphere`.
"""
struct raMOInput
    dftdata::raMODFTData
    runlist::Vector{RunInfo}
    emin::Float64
    emax::Float64
    checkpoint::String  # TODO: should we use some sort of IO type? or checkpoint container?
    auto_psphere::Bool
    discard::Bool
    function raMOInput(
        dftdata::raMODFTData,
        runlist,
        emin::Real,
        emax::Real,
        checkpoint::AbstractString,
        auto_psphere,
        discard
    )
        @assert emin < emax "Minimum energy is not less than maximum energy"
        isnothing(auto_psphere) && (auto_psphere = false)
        isnothing(discard) && (discard = false)
        return new(dftdata, collect(runlist), emin, emax, checkpoint, auto_psphere, discard)
    end
end

"""
    raMOInput(
        dftdata::raMODFTData,
        runlist;
        auto_psphere = false,
        checkpoint = "",
        emin = minimum(dftdata.wave.energies),
        emax = fermi(dftdata)
    )

Constructs a `raMOInput` object with some assumptions implemented as keyword arguments.

`auto_psphere` is set to false if not specified.

If `checkpoint` is unset, the checkpoint path is the empty string, corresponding to no checkpoint
file being used.

If `emin` is unset, then the value is automatically set to the lowest energy in the range of the
`PlanewaveWavefunction` within the `raMODFTData`.

If `emax` is unset, then the value is automatically set to the Fermi energy of the wavefunction.
"""
function raMOInput(
    dftdata::raMODFTData,
    runlist;
    auto_psphere = false,
    discard = false,
    checkpoint::AbstractString = "",
    emin::Real = minimum(dftdata.wave.energies),
    emax::Real = fermi(dftdata)
)
    return raMOInput(dftdata, runlist, emin, emax, checkpoint, auto_psphere, discard)
end

raMODFTData(x::raMOInput) = x.dftdata

Electrum.Crystal(x::raMOInput) = x.dftdata.xtal
Electrum.basis(x::raMOInput) = basis(x.dftdata.xtal)
Electrum.PlanewaveWavefunction(x::raMOInput) = x.dftdata.wave
Electrum.PeriodicAtomList(x::raMOInput) = x.dftdata.xtal.atoms
Electrum.supercell(x::raMOInput) = supercell(x.dftdata)

# Size and indexing depend on the run list.
Base.size(x::raMOInput) = size(x.runlist)
Base.getindex(x::raMOInput, i) = x.runlist[i]

"""
    Supercell(atomlist::Atomlist{3}, orblist_by_type::Dict{String,Int64})

Returns a supercell struct AtomList{3} and a Vector{Int} with number of orbitals associated with each atom
"""
struct Supercell
    atomlist::PeriodicAtomList{3}
    orbitals::Vector{Int}
    function Supercell(ramoinput::raMOInput, orblist_by_type::Dict{String,Int64})
        orbitals = Vector{Int}(undef,0)
        atomlist = supercell(ramoinput)
        for atom in atomlist
            if !haskey(orblist_by_type, name(atom))
                error("Number of orbitals for atom", name(atom),"cannot be found.")
            end
            # Get the number of orbitals for current atom and add to sum
            push!(orbitals,get(orblist_by_type, name(atom),0))
        end
        return new(atomlist, orbitals)
    end
end

Electrum.basis(s::Supercell) = basis(s.atomlist)
Electrum.PeriodicAtomList(s::Supercell) = s.atomlist
Base.getindex(s::Supercell, i...) = getindex(s.atomlist, i...)
