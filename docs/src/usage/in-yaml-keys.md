# Keys for the Instructional YAML

## General settings in the instructional YAML file
| Key | Purpose | Options |
| :--- | :--- | :--- |
| `checkpoint` | Allows the current sequence of DFT-raMO to begin using the specified basis set typically obtained partway through the analysis. | Blank – starts with the basis set generated from DFT inputs |
| | | Path to the .chkpt file containing the starting basis set for the analysis |
| `mode` | Specifies if DFT-raMO evaluates the raMO and how the basis set is updated. | `default` – keeps all generated raMOs regardless of ``P_{sphere}`` value and updates the basis set accordingly. |
| | | `auto_psphere` – automatically returns raMOs with ``P_{sphere}`` values <15% of the maximum ``P_{sphere}`` value to the basis set. Subsequent raMOs that do meet the criteria are automatically rerun. |
| | | `discard` – all raMOs are returned to the basis set. Typically used for testing raMOs at different sites. |
| `software` | Indicates which DFT package is used to calculate wavefunction information | `vasp` (`abinit` is currently not supported but will be implemented in a future update) |
| `emin` and `emax` | Specifies the lower and upper limit of energy for wavefunctions used in the basis set | Blank – by default, will start at the lowest energy |
| | |	(Number) eV – Energy in eV |
| | | (Number) Ha – Energy in Ha |
| `runs` | Indicates the list of raMO sequences will follow | The list of raMO sequences and corresponding keys. See [table](#*runs*-list). |

## *runs* list
| Key | Purpose | Options |
| :--- | :--- | :--- |
| `name` | Names the directory that stores output files corresponding to the raMO analysis. Also prepends files within that directory with this name. |	Blank – (not recommended) by default, starts with `run_1` and continues numbering. |
| | | A string. It is recommended to use a descriptive yet concise string with no whitespace (e.g. `1_Sc_dxy` to indicate it is the first sequence and corresponds to reconstruction of the Sc dxy atomic orbitals). |
| `type` | Specifies which type of target raMO to generate | Atomic orbitals (`s`, `px`, `py`, `pz`, `dx2y2` or `dx2-y2`, `dz2`, `dxy`, `dxz`, `dyz`) |
| | | `sp` - Molecular orbital-like targets built from *s*- and *p*-based orbitals based on distance to a specified central site. |
| | | `lcao` - Builds targets based on user-specified linear combination of atomic orbitals. |
| | | `displaced (AO)` - Creates targets that are displaced off the central site by some specified direction and distance. |
| `site_file` | Specifies the .xyz (*sp* type) or .yaml (*lcao* type) file for the sequence | For *sp* type runs, the path to an .xyz file |
| | | For *lcao* type runs, the path to a .yaml file |
| `sites` | Specifies which sites in the periodic atom list/.xyz/.yaml to build targets for | For atomic orbital type runs, a list of atomic sites in the supercell based on the periodic atom list (e.g. `1, 2, 4:12`) or specific elements (e.g. `Sc` for all Sc sites, or `Sc 1, 2, 4:12` for the 1st, 2nd, and 4th through 12th sites of Sc) |
| | | For *sp* and *lcao* type runs, `all` builds targets for all sites listed in the .xyz or .yaml file, or a list of sites based on the .xyz or .yaml file (e.g. `1, 2, 4:12`) |
| `direction` | For *displaced (AO)* type runs only, specifies a directional vector in fractional coordinates for the displacement | `[a, b, c]` – where a, b, and c correspond to the displacement direction in the fractional coordinate system. The directional vector is automatically normalized. |
| `radius` | For *sp* type runs, indicates the maximum distance from the site to search for atoms whose *s* and *p* orbitals are used in generating the target. For *displaced (AO)* type runs, indicates the distance of the displacement. | `(Number)` - distance in Å |
| `rsphere` | Specifies ``r_{sphere}`` used for ``P_{sphere}`` analysis | Blank - a default value of 3 Å will be applied |
| | | `(Number)` - distance in Å |