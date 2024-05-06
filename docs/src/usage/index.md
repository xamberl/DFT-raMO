# Overview

DFT-raMO analysis occurs in multiple steps. First, a DFT package is used to generate the crystal orbitals of the system. Then DFT-raMO.jl imports the information and reads a YAML file that contains user-defined instructions for the analysis. DFT-raMO.jl will then execute the analysis sequentially, using the crystal orbitals as a basis set to generate MO-like functions (raMOs). These raMOs can be quickly assessed through the ``P_{sphere}`` metric and inform the user on the validity and viability of the proposed bonding states.

Jump to:
- [DFT calculations](DFT-calculations.md)
- [The Instructional YAML file](in-yaml.md)
  - [LCAO sites file](LCAO.md)
  - [Keys for the Instructional YAML](in-yaml-keys.md)

# Running DFT-raMO
After configuring the basic input file, ensure that you are in your working directory, which includes
the `OUTCAR`, `POSCAR`, `KPOINTS`, and `WAVECAR` files. Open the Julia REPL.
```julia
using DFTraMO
dftramo_run("input.yaml")
```
