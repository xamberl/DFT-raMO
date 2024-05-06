# Example system: IrIn3

In this example, we're going to demonstrate a DFT-raMO run for an intermetallic binary compound in
the iridium-indium system, using data from a VASP 4.6 calculation.

## Setup

Assuming that you have Julia and DFTraMO.jl installed, download the IrIn3 VASP 4.6 data, a link to
which will be provided in the future. This includes the following files:
  
  - `KPOINTS`
  - `OUTCAR`
  - `POSCAR`
  - `POTCAR`
  - `WAVECAR`
  - `step1.yaml` - the DFT-raMO input file.
  - `qdftramo` and `qdftramo-kestrel` - submission scripts for internal use, can be adapted.

!!! warn The tarball is a large file due to the `WAVECAR` (over 400 MB).

## Constructing the YAML input

You'll find that `step1.yaml` is a blank template. The [Usage] section of the manual contains
information on each of the keys and values in the file.

Since we're starting from scratch, the `checkpoint` key should be blank. Leave `auto_psphere` as
`true`, as this will automatically rerun poor reconstructions (less than 15% of the maximum
``P_sphere`` value).

For the first run, let's reconstruct the Ir ``d_{x^2+y^2}`` orbitals, so set `type` to `dx2y2` and
`sites` to `Ir`. `site_file` and `radius` can be blank, as these are only used for AO type runs.
As a rule of thumb, you can use the distance from the reconstruction site to the farthest atom in
the first coordination sphere, which is 2.65 angstroms (this does not have to be precise, but
overshooting is generally better than undershooting).

The `name` field can be set at your convenience (here we use `1_Ir_dx2y2`) and becomes the name of
the directory where the outputs are written.

## Running DFT-raMO for the first time

In the directory you extracted the tarball, run julia, use the DFTraMO.jl package, and start the
run:

```julia-repl
julia> using DFTraMO

julia> dftramo_run("step1.yaml")
```

Optionally, you can use the `UnicodePlots` package to visualize the Psphere data. If you don't load
the package

## Examining the output

## Reconstructing sp hybrids

In order to reconstruct sp hybrids, we want to reconstruct all of the Ir d and s orbitals. Instead
of running all of the reconstructions yourself, to save time, the `step2.tar` file contains all of
the necessary reconstructions.

For this next run, create a new YAML input with the following lines:
```yaml
checkpoint: 6_Ir_s/6_Ir_s_192_192.chkpt
auto_psphere: true
runs:
  - name: 7_Ir-Ir
    type: sp
    site_file: Ir-Ir.txt
    sites: all
    radius: 1.6
    rsphere: 3.5
```
Our previous analysis ended at the 192nd electron, and we have 192 electrons remaining, explaining
the occurrence of 192 twice in the name of the checkpoint file. The `radius` parameter is the search
radius for automatically including atoms in the target. 

The `Ir-Ir.txt` site file should contain the midpoints of the Ir-Ir bonds in Cartesian coordinates.
This can be done by using your preferred molecular modeling software (we've done this by converting
the POSCAR to a CIF and adding the dummy atoms into Diamond 3). However, if you'd like to skip the
manual creation, the contents of `Ir-Ir.txt` are given below:
```
X	  6.99330	  6.99330	 10.78620
X	  0.00000	  6.99330	 10.78620
X	  6.99330	  0.00000	 10.78620
X	  0.00000	  0.00000	 10.78620
X	  6.99330	  6.99330	  3.59540
X	  6.99330	  0.00000	  3.59540
X	  0.00000	  0.00000	  3.59540
X	 10.48995	 10.48995	  7.19080
X	  3.49665	 10.48995	  7.19080
X	 10.48995	  3.49665	  7.19080
X	  0.00000	  6.99330	  3.59540
X	  3.49665	  3.49665	  7.19080
X	 10.48995	 10.48995	  0.00000
X	  3.49665	 10.48995	  0.00000
X	 10.48995	  3.49665	  0.00000
X	  3.49665	  3.49665	  0.00000
```
The symbol `X` can be changed arbitrarily, as long as it does not match the symbol of another atom
already present.

Now that you have these files, you can just run DFT-raMO again:
```julia-repl
julia> dftramo_run("step2.yaml")
```

As before, have a look at the XSF files to see the reconstructions.

## LCAO reconstructions

When working with orbitals with nonzero orbital angular momentum (p, d, f orbitals) it may be
desirable to construct an orbital with an orientation or shape that is not aligned with the
coordinate system. To solve this issue, we can use linear combinations of atomic orbitals (LCAOs).

To skip to the next step, untar `step3.tar` in your working directory. You should find two new input
files, `lcao1.yaml` and `lcao2.yaml`. The contents of `lcao1.yaml` are given below:
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
As before, the description of the LCAO sites file is given in the [Usage] section.

The `lcao2.yaml` file contains reconstructions with different orientations

The `step3.yaml` file contains the information for the LCAO reconstruction, and references both
`lcao1.yaml` and `lcao2.yaml` for each run. Its contents are below:
```yaml
checkpoint: 8_Ir_pz/8_Ir_pz_240_96.chkpt
auto_psphere: true
runs:
  - name: 9_Ir_p_pi1
    type: lcao
    site_file: lcao1.yaml
    sites: all
    radius:
    rsphere: 2.65
  - name: 10_Ir_p_pi2
    type: lcao
    site_file: lcao2.yaml
    sites: all
    radius:
    rsphere: 2.65
```

# Remainder analysis

As you complete your DFT-raMO runs, you'll find that the Psphere values of the reconstructed
orbitals tends to decrease. To visualize this, untar `step4.tar` and look at the log for run 11,
found at `11_In-In/11_In-In_psphere_4.73.txt`:
```
273     0.6474487830205945       at site [10.490, 3.497, 3.595]
274     0.6429286864719969       at site [0.000, 6.993, 7.191]
275     0.6116188776692522       at site [3.497, 10.490, 3.595]
276     0.6382205090625862       at site [10.490, 10.490, 3.595]
277     0.6504732377051952       at site [6.993, 6.993, 7.191]
278     0.5981010110158822       at site [3.497, 3.497, 10.786]
279     0.5871384506957702       at site [10.490, 3.497, 10.786]
280     0.5162750706706059       at site [3.497, 10.490, 10.786]
281     0.33634765527327326      at site [10.490, 10.490, 10.786]
282     0.6434855799801233       at site [0.000, 0.000, 0.000]
283     0.6433852489950562       at site [6.993, 0.000, 0.000]
284     0.5527474698233648       at site [0.000, 6.993, 0.000]
285     0.49587592534638286      at site [6.993, 6.993, 0.000]
286     0.40417825671730373      at site [0.000, 0.000, 7.191]
287     0.16723350141820606      at site [6.993, 0.000, 7.191]
288     0.1572312564858335       at site [3.497, 3.497, 3.595]
```
We find that ``P_{sphere}`` drops precipitously at some points in this run. In order to complete the
analysis, we want to remove the orbitals from our total reconstruction and 
We recommend that you manually inspect the orbitals to observe their
character. When we do this, we find that orbitals 287 and 288 do not have the desired character, and
the last run should reconstruct from orbital 286.

For the final run, we can create the last input file, `step4.yaml`:
```yaml
checkpoint: 11_In-In/11_In-In_286_4.chkpt
auto_psphere: true
runs:
  - name: 12.1_In_px
    type: px
    site_file: 
    sites: In 1
    radius:
    rsphere: 3.2
  - name: 12.2_In_py
    type: py
    site_file: 
    sites: In 1
    radius:
    rsphere: 3.2
  - name: 12.3_In_pz
    type: pz
    site_file: 
    sites: In 1
    radius:
    rsphere: 3.2
```
