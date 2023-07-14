"""
    raMOInput

Contains all of the crystal and wavefunction information needed to perform a DFT-raMO run.
"""
struct raMOInput
    xtal::Crystal
    wave::PlanewaveWavefunction
    fermi::Float64
end

Electrum.basis(x::raMOInput) = basis(x.xtal.atoms)
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

"""
    raMOInput(POSCAR, WAVECAR, KPOINTS, OUTCAR, ::FromVASP)
    raMOInput(;
        POSCAR = "POSCAR",
        WAVECAR = "WAVECAR",
        KPOINTS = "KPOINTS",
        OUTCAR = "OUTCAR"
    )

Reads VASP output files required for DFT-raMO. The files may be specified as strings, I/O handles,
or path types if another package provides them (for more details, see the Electrum functions
`readPOSCAR`, `readWAVECAR`, `readKPOINTS`, and `get_fermi`).

File names may be specified through the argument order given above, or with keyword arguments which
have no fixed order.
"""
function raMOInput(POSCAR, WAVECAR, KPOINTS, OUTCAR, ::FromVASP)
    fermi = get_fermi(OUTCAR) * Electrum.EV2HARTREE
    geo = readPOSCAR(POSCAR)
    wave = readWAVECAR(WAVECAR, quiet = true)
    kpt = parse.(Int, split(readlines(KPOINTS)[4]))
    # Use a Crystal to lazily reference the supercell
    xtal = set_transform!(Crystal(geo), kpt)
    return raMOInput(xtal, wave, fermi) 
end

function raMOInput(;
    POSCAR = "POSCAR",
    WAVECAR = "WAVECAR",
    KPOINTS = "KPOINTS",
    OUTCAR = "OUTCAR"
)
    return raMOInput(POSCAR, WAVECAR, KPOINTS, OUTCAR, FromVASP())
end

"""
    raMOInput(directory, FromVASP(); CONTCAR=false)
    raMOInput(FromVASP(); CONTCAR=false)

Reads VASP output files required for DFT-raMO from the given directory. If no directory is given,
the current directory is checked.
    
If `CONTCAR` is set to true, the CONTCAR file will be used instead of the POSCAR.
"""
function raMOInput(directory, ::FromVASP; CONTCAR = false)
    return raMOInput(
        joinpath(directory, CONTCAR ? "CONTCAR" : "POSCAR"),
        joinpath(directory, "WAVECAR"),
        joinpath(directory, "KPOINTS"),
        joinpath(directory, "OUTCAR"),
        FromVASP()
    )
end

raMOInput(::FromVASP; CONTCAR) = raMOInput(".", FromVASP(); CONTCAR)
