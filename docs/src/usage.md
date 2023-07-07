# Usage

## Input files

Currently, DFT-raMO supports VASP inputs, and the required VASP outputs are the `OUTCAR`, `POSCAR`,
`KPOINTS`, and `WAVECAR` files. Future versions will support abinit (likely to be tested on versions
8.10.3 and 9.10.1), and this will likely require a `WFK` output as well as the NetCDF `OUT.nc` file.

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
