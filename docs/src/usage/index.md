# Input files

## Preparation of VASP input files
Currently, DFT-raMO supports VASP inputs, and the required VASP outputs are the `OUTCAR`, `POSCAR`,
`KPOINTS`, and `WAVECAR` files. Future versions will support abinit (likely to be tested on versions
8.10.3 and 9.10.1), and this will require a `WFK` output.

A static calculation should be performed with a `POSCAR` pertaining to the unit cell geometry. The `POSCAR` file must specify atomic identities after giving cell parameters but before the stoichiometry and direct coordinates. This is not a problem with VASP version 5 or later, but POSCAR files in VASP 4 do not include this information. Ensure this information is in the `POSCAR` or add it manually. An example is given below, where the atomic identities are given in the line `Sc Al`.
```
Al3 Sc1                                 
   1.00000000000000     
     4.1046969868448775    0.0000000000000000    0.0000000000000000
     0.0000000000000000    4.1046969868448775    0.0000000000000000
     0.0000000000000000    0.0000000000000000    4.1046969868448775
Sc Al
   1   3
Direct
  0.0000000000000000  0.0000000000000000  0.0000000000000000
  0.0000000000000000  0.5000000000000000  0.5000000000000000
```

A Γ-centered k-point mesh (`KPOINTS`) corresponding to the desired supercell in real space should also be used. For example, if using a 2×2×3 supercell, the k-point grid must also be 2×2×3. The supercell size should be large enough such that atoms do not interact with their translational copies in neighboring supercells. The recommended starting point is to use a supercell approximately 12 Å in length for each axis (assuming cubic, tetragonal, or orthorhombic systems). A larger grid may be used but the user should be aware that significantly more memory will be needed.

Be sure to turn off symmetry (`ISYM = -1`) and print the `WAVECAR`.


## The Instructional YAML file

!!! note
    See also [Keys for the Instructional YAML](in-yaml-keys.md)

DFT raMO.jl requires an instructions file to customize and perform the analysis. This file specifies input files, checkpoint files, and other necessary parameters. The input file must be written in YAML format (with the .yml or .yaml extension), which consists of a series of keys and values separated by colons. Lists are defined with a hyphen and spaces. Below is an example of inputs for the instructions file for DFT-raMO.jl.

Example YAML file:
```yaml
checkpoint:
software: vasp
mode: auto_psphere
emin: -100 eV
emax:
runs:
   - name: 1_Ir_dx2y2
     type: dx2y2
     sites_file: 
     sites: Ir
     rsphere: 2.65
   - name: 2_Ir_dz2
     type: dz2
     sites_file: 
     sites: Ir
     rsphere: 2.65
```

Here are some of the important components of this example file:
  - The `checkpoint` key can be left blank when starting a new calculation to use an initial basis set generated from the DFT calculations. Checkpoint files are automatically written by DFT-raMO. A path to a file can be given to resume the calculation. 
  - The `software` key specifies which DFT software generated the wavefunction information.
  - The `mode` key is set to `auto_psphere` mode, allowing DFT-raMO to automatically reject functions that do not meet the criteria set by ``P_{sphere}`` analysis.
  - The `emin` and `emax` keys allow the user to specify the energy range (in eV or Ha) of bands used in the initial basis set — leaving these values blank will default to using all bands below the Fermi energy.
  - The `runs` key indicates a list of raMO sequences will follow. Each sequence is prepended with spaces and a hyphen. Keys within the list item are exclusive to that sequence.
      + The `name` key requires a string which will be used to ame the directory that stores the sequence's output files. By default, the name is `run_<number>` where `<number>` is the number in the order.
      + The `type` key determines the target orbital shapes to reconstruct: atomic orbitals (`s`, `px`, `py`, `pz`, `dx2y2` or `dx2-y2`, `dz2`, `dxy`, `dxz`, and `dyz`), sp-based orbitals built from a distance criteria (`sp`), linear combination of atomic orbitals (`lcao`), or displaced atomic orbitals (`displaced (AO)`).
      + The `rsphere` key corresponds to the distance from the central site to consider for ``P_{sphere}`` analysis in Angstroms.

A complete table of general keys and options are listed in the [Keys for the Instructional YAML](in-yaml-keys.md).

## LCAO sites file
Here is an example YAML file for user-defined linear combinations of atomic orbitals (LCAO).

```yaml
target:
  - px: -1
    py: 1
lcao:
  - [1]
  - [2]
  - [17]
  - [18]
  - [33]
  - [34]
  - [49]
  - [50]
  - [65]
  - [66]
  - [81]
  - [82]
  - [97]
  - [98]
  - [113]
  - [114]
```

Here are some of the important components of this example file:
  - The `target:` key indicates that a list of LCAOs will follow. Each atomic site is indicated with
    a `-`, and relative contributions of each atomic orbital follow each atomic orbital key.
      + Relative contributions of each atomic orbital will be normalized such that the total
        contribution equals 1.
      + Any atomic orbital keys that are notdefined are ignored.
  - The `lcao` key indicates the list of targets in the run.
      + Each list item must be a vector of integers indicating the number corresponding to the atom
        in the supercell. These atomic positions can usually be found by checking the xsfs of a
        previously run.
      + The number of integers in the list *must* match the number of atomic sites in `target`.

# Running DFT-raMO
After configuring the basic input file, ensure that you are in your working directory, which includes
the `OUTCAR`, `POSCAR`, `KPOINTS`, and `WAVECAR` files. Open the Julia REPL.
```julia
using DFTraMO
dftramo_run("input.yaml")
```
