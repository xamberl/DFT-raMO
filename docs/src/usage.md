# Input files

## Basic input file

Currently, DFT-raMO supports VASP inputs, and the required VASP outputs are the `OUTCAR`, `POSCAR`,
`KPOINTS`, and `WAVECAR` files. Future versions will support abinit (likely to be tested on versions
8.10.3 and 9.10.1), and this will require a `WFK` output.

In most cases, you will want to use an optimized geometry and a gamma-centered k-point mesh which
corresponds to the size of the supercell you want to reconstruct the orbitals in.

Regardless of the software used to generate the wavefunction, DFT-raMO uses a YAML file as input
which contains information about the desired target reconstruction.

Here is an example YAML file (available at `examples/ScAl3.yaml` and `test/input/ScAl3.yaml`):

```yaml
checkpoint:
auto_psphere: true
runs:
  - name: 1_Sc-Sc
    type: sp
    site_file: Sc-Sc_reordered.xyz
    sites: all
    radius: 2.2
    rsphere: 2.2
  - name: 2_Sc_dxy
    type: dxy
    site_file: 
    sites: Sc
    radius: 
    rsphere: 2.2
  - name: 3_Sc_dxz
    type: dxz
    site_file: 
    sites: Sc
    radius: 
    rsphere: 2.2
  - name: 4_Sc_dyz
    type: dyz
    site_file: 
    sites: Sc
    radius: 
    rsphere: 2.2
  - name: 5_Al_oct
    type: sp
    site_file: Al-Al_reordered.xyz
    sites: all
    radius: 2.1
    rsphere: 4.11
```

Here are some of the important components of this example file:
  - The `checkpoint` key can be left blank when starting a new calculation, but a path to a file can
  be given to resume the calculation. Checkpoint files are automatically written.
  - The `auto_psphere` key automatically reruns calculations where the Psphere values are less than
    15% of the maxium when set to `true`.
  - The `runs` key is the bread and butter of the sequence.
      + The `name` key can be used to name the output directory and the output files it creates. By
        default, the name is `run_<number>` where `<number>` is the number in the order.
      + The `type` key determines the target orbital shapes to reconstruct. By default, it accepts
        atomic orbital designations (`s`, `px`, `py`, `pz`, `dx2y2` or `dx2-y2`, `dz2`, `dxy`, 
        `dxz`, and `dyz`), or hybrid designations ( `sp`, which includes p-orbitals
        with any directionality, or user-defined linear combinations of atomic orbitals with `lcao`)
        . These designations are not case-sensitive.
      + `site_file` is only used for `sp` cage and `lcao` runs. The input is an
        XYZ file (Cartesian coordinates in angstroms).
      + `sites` can accept multiple arguments. By default, this is `all`, which runs the whole list
        of atoms. However, this can be limited if desired: the atom names from the POSCAR can be
        included, and can be postpended with lists of individual sites or ranges of sites.
      + The `radius` is only relevant for `sp` hybrids, and it is the length (in angstroms) used to
        search for atoms that contribute to that hybrid.
      + `rsphere` is the radius for the Psphere metric.


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
