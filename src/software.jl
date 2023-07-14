"""
    raMOInput

Contains all of the crystal and wavefunction information needed to perform a DFT-raMO run.
"""
struct raMOInput
    xtal::Crystal
    wave::PlanewaveWavefunction
    fermi::Float64
end

Electrum.Crystal(x::raMOInput) = x.xtal
Electrum.PeriodicAtomList(x::raMOInput) = x.xtal.atoms
Electrum.PlanewaveWavefunction(x::raMOInput) = x.wave
Electrum.fermi(x::raMOInput) = x.fermi
kptmesh(x::raMOInput) = x.xtal.set_transform

"""
    DFTraMO.InputOrigin

Dispatch type for various computational chemistry packages (for a complete list, run
`subtypes(DFTraMO.InputOrigin)` in the REPL).
"""
abstract type InputOrigin
end

"""
    DFTraMO.FromABINIT

Dispatch type for reading abinit WFK outputs.
"""
struct FromABINIT
end

"""
    DFTraMO.FromVASP

Dispatch type for reading VASP calculation outputs (specifically, the `POSCAR`, `WAVECAR`,
`KPOINTS`, and `OUTCAR` files.)
"""
struct FromVASP
end

function raMOInput(io::IO, ::FromABINIT)
    h = Electrum.read_abinit_header(io)
    seekstart(io)
    return raMOinput(Crystal(h), read_abinit_WFK(io)["wavefunction"], h[:fermi])
end

raMOInput(file, ::FromABINIT) = open(io -> raMOInput(io, FromABINIT()), file)

function raMOInput(
    ::FromVASP;
    POSCAR = "POSCAR",
    WAVECAR = "WAVECAR",
    KPOINTS = "KPOINTS",
    OUTCAR = "OUTCAR"
)
    fermi = get_fermi(OUTCAR)
    geo = readPOSCAR(POSCAR)
    wave = readWAVECAR(WAVECAR, quiet = true)
    kpt = parse.(Int, split(readlines(KPOINTS)[4]))
    # Use a Crystal to lazily reference the supercell
    xtal = set_transform!(Crystal(geo), kpt)
    return raMOInput(xtal, wave, fermi) 
end

function raMOInput(directory::AbstractString, ::FromVASP)
    return raMOInput(
        FromVASP();
        POSCAR = joinpath(directory, POSCAR),
        WAVECAR = joinpath(directory, WAVECAR),
        KPOINTS = joinpath(directory, KPOINTS),
        OUTCAR = joinpath(directory, OUTCAR),
    )
end
