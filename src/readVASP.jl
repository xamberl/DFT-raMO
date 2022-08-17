# Initial script to read VASP inputs
using Xtal

# Transform WAVECAR to data [TBD!!!]

# Grab fermi energy from OUTCAR
function getFermi(io::IO)
    readuntil(io, "E-fermi :")
    fermi = parse.(Float64, split(readline(io))[1])
    return fermi
end

# Grab supercell dimensions from KPOINTS
function getSupercell(io::IO)
    supercelldims = parse.(Int64,split(readlines(io)[4]))
    return supercelldims
end

export getFermi, getSupercell

