"""
    DFTraMO.InputOrigin{S}

Dispatch type to indicate the software package which generated the input files for DFT-raMO. The
type parameter `S` is a `Symbol`, all lowercase, which contains the name of the software package.

To see a list of the supported software packages, evaluate `DFTraMO.SUPPORTED_SOFTWARE`.
"""
struct InputOrigin{S}
end

const SUPPORTED_SOFTWARE = (:abinit, :vasp)

"""
    DFTraMO.length_conversion(i::InputOrigin) -> Real

Provides the conversion factor for the default length units for the software package indicated by
`i` to Bohr. By default, this is equal to 1 (we assume most packages use Hartree atomic units).
"""
length_conversion(::InputOrigin) = 1

"""
    DFTraMO.energy_conversion(i::InputOrigin) -> Real

Provides the conversion factor for the default energy units for the software package indicated by
`i` to Hartree. By default, this is equal to 1 (we assume most packages use Hartree atomic units).
"""
energy_conversion(::InputOrigin) = 1

"""
    DFTraMO.FromABINIT

Dispatch type for reading abinit WFK outputs.
"""
const FromABINIT = InputOrigin{:abinit}

"""
    DFTraMO.FromVASP
 
Dispatch type for reading VASP calculation outputs (specifically, the `POSCAR`, `WAVECAR`,
`KPOINTS`, and `OUTCAR` files.)
"""
const FromVASP = InputOrigin{:vasp}
# Default units for VASP are angstroms and electron-volts
energy_conversion(::FromVASP) = Electrum.EV2HARTREE
length_conversion(::FromABINIT) = Electrum.BOHR2ANG

function raMODFTData(io::IO, ::FromABINIT)
    h = Electrum.read_abinit_header(io)
    seekstart(io)
    return raMODFTData(Crystal(h), read_abinit_WFK(io)["wavefunction"], h[:fermi])
end

raMODFTData(file, ::FromABINIT) = open(io -> raMODFTData(io, FromABINIT()), file)

"""
    raMODFTData(POSCAR, WAVECAR, KPOINTS, OUTCAR, ::FromVASP)
    raMODFTData(;
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
function raMODFTData(POSCAR, WAVECAR, KPOINTS, OUTCAR, ::FromVASP)
    fermi = get_fermi(OUTCAR).fermi
    geo = readPOSCAR(POSCAR)
    wave = readWAVECAR(WAVECAR, quiet = true)
    kpt = parse.(Int, split(readlines(KPOINTS)[4]))
    # Use a Crystal to lazily reference the supercell
    xtal = set_transform!(Crystal(geo), kpt)
    return raMODFTData(xtal, wave, fermi) 
end

function raMODFTData(;
    POSCAR = "POSCAR",
    WAVECAR = "WAVECAR",
    KPOINTS = "KPOINTS",
    OUTCAR = "OUTCAR"
)
    return raMODFTData(POSCAR, WAVECAR, KPOINTS, OUTCAR, FromVASP())
end

"""
    raMODFTData(directory, FromVASP(); CONTCAR=false)
    raMODFTData(FromVASP(); CONTCAR=false)

Reads VASP output files required for DFT-raMO from the given directory. If no directory is given,
the current directory is checked.
    
If `CONTCAR` is set to true, the CONTCAR file will be used instead of the POSCAR.
"""
function raMODFTData(directory, ::FromVASP; CONTCAR = false)
    return raMODFTData(
        joinpath(directory, CONTCAR ? "CONTCAR" : "POSCAR"),
        joinpath(directory, "WAVECAR"),
        joinpath(directory, "KPOINTS"),
        joinpath(directory, "OUTCAR"),
        FromVASP()
    )
end

raMODFTData(::FromVASP; CONTCAR) = raMODFTData(".", FromVASP(); CONTCAR)
