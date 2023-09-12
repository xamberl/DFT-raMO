var documenterSearchIndex = {"docs":
[{"location":"usage.html#Input-files","page":"Usage","title":"Input files","text":"","category":"section"},{"location":"usage.html#Basic-input-file","page":"Usage","title":"Basic input file","text":"","category":"section"},{"location":"usage.html","page":"Usage","title":"Usage","text":"Currently, DFT-raMO supports VASP inputs, and the required VASP outputs are the OUTCAR, POSCAR, KPOINTS, and WAVECAR files. Future versions will support abinit (likely to be tested on versions 8.10.3 and 9.10.1), and this will require a WFK output.","category":"page"},{"location":"usage.html","page":"Usage","title":"Usage","text":"Regardless of the software used to generate the wavefunction, DFT-raMO uses a YAML file as input which contains information about the desired target reconstruction.","category":"page"},{"location":"usage.html","page":"Usage","title":"Usage","text":"Here is an example YAML file (available at examples/ScAl3.yaml and test/input/ScAl3.yaml):","category":"page"},{"location":"usage.html","page":"Usage","title":"Usage","text":"checkpoint:\nauto_psphere: true\nruns:\n  - name: 1_Sc-Sc\n    type: sp\n    site_file: Sc-Sc_reordered.xyz\n    sites: all\n    radius: 2.2\n    rsphere: 2.2\n  - name: 2_Sc_dxy\n    type: dxy\n    site_file: \n    sites: Sc\n    radius: \n    rsphere: 2.2\n  - name: 3_Sc_dxz\n    type: dxz\n    site_file: \n    sites: Sc\n    radius: \n    rsphere: 2.2\n  - name: 4_Sc_dyz\n    type: dyz\n    site_file: \n    sites: Sc\n    radius: \n    rsphere: 2.2\n  - name: 5_Al_oct\n    type: sp\n    site_file: Al-Al_reordered.xyz\n    sites: all\n    radius: 2.1\n    rsphere: 4.11","category":"page"},{"location":"usage.html","page":"Usage","title":"Usage","text":"Here are some of the important components of this example file:","category":"page"},{"location":"usage.html","page":"Usage","title":"Usage","text":"The checkpoint key can be left blank when starting a new calculation, but a path to a file can","category":"page"},{"location":"usage.html","page":"Usage","title":"Usage","text":"be given to resume the calculation. Checkpoint files are automatically written.","category":"page"},{"location":"usage.html","page":"Usage","title":"Usage","text":"The auto_psphere key automatically reruns calculations where the Psphere values are less than 15% of the maxium when set to true.\nThe runs key is the bread and butter of the sequence.\nThe name key can be used to name the output directory and the output files it creates. By default, the name is run_<number> where <number> is the number in the order.\nThe type key determines the target orbital shapes to reconstruct. By default, it accepts atomic orbital designations (s, px, py, pz, dx2y2 or dx2-y2, dz2, dxy,  dxz, and dyz), or hybrid designations ( sp, which includes p-orbitals with any directionality, or user-defined linear combinations of atomic orbitals with lcao) . These designations are not case-sensitive.\nsite_file is only used for sp cage and lcao runs. The input is an XYZ file (Cartesian coordinates in angstroms).\nsites can accept multiple arguments. By default, this is all, which runs the whole list of atoms. However, this can be limited if desired: the atom names from the POSCAR can be included, and can be postpended with lists of individual sites or ranges of sites.\nThe radius is only relevant for sp hybrids, and it is the length (in angstroms) used to search for atoms that contribute to that hybrid.\nrsphere is the radius for the Psphere metric.","category":"page"},{"location":"usage.html#LCAO-sites-file","page":"Usage","title":"LCAO sites file","text":"","category":"section"},{"location":"usage.html","page":"Usage","title":"Usage","text":"Here is an example YAML file for user-defined linear combinations of atomic orbitals (LCAO).","category":"page"},{"location":"usage.html","page":"Usage","title":"Usage","text":"target:\n  - px: -1\n    py: 1\nlcao:\n  - [1]\n  - [2]\n  - [17]\n  - [18]\n  - [33]\n  - [34]\n  - [49]\n  - [50]\n  - [65]\n  - [66]\n  - [81]\n  - [82]\n  - [97]\n  - [98]\n  - [113]\n  - [114]","category":"page"},{"location":"usage.html","page":"Usage","title":"Usage","text":"Here are some of the important components of this example file:","category":"page"},{"location":"usage.html","page":"Usage","title":"Usage","text":"The target: key indicates that a list of LCAOs will follow. Each atomic site is indicated with a -, and relative contributions of each atomic orbital follow each atomic orbital key.\nRelative contributions of each atomic orbital will be normalized such that the total contribution equals 1.\nAny atomic orbital keys that are notdefined are ignored.\nThe lcao key indicates the list of targets in the run.\nEach list item must be a vector of integers indicating the number corresponding to the atom in the supercell. These atomic positions can usually be found by checking the xsfs of a previously run.\nThe number of integers in the list must match the number of atomic sites in target.","category":"page"},{"location":"usage.html#Running-DFT-raMO","page":"Usage","title":"Running DFT-raMO","text":"","category":"section"},{"location":"usage.html","page":"Usage","title":"Usage","text":"After configuring the basic input file, ensure that you are in your working directory, which includes the OUTCAR, POSCAR, KPOINTS, and WAVECAR files. Open the Julia REPL. julia using DFTraMO dftramo_run(\"input.yaml\")```","category":"page"},{"location":"theory.html#Theory","page":"Theory","title":"Theory","text":"","category":"section"},{"location":"theory.html","page":"Theory","title":"Theory","text":"Coming soon...","category":"page"},{"location":"api/index.html#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/index.html","page":"API","title":"API","text":"Coming soon...","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"(Image: DFT-raMO.jl logo)","category":"page"},{"location":"index.html#DFT-raMO.jl","page":"Home","title":"DFT-raMO.jl","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"The DFT-raMO (reversed approximation Molecular Orbital) method is a Wannier-type analysis that can be used to reinterpret the results of an electronic structure calculation in terms of localized orbitals, similar to natural hybrid orbital (NHO) and natural bond orbital (NBO) analysis.","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"Unlike NBO, raMO follows user-guided targets, which can be centered not just on atoms, but on any user-selected site. The reconstructions are not dependent on the basis set used to perform the calculation, which makes this method ideal for solid-state structures.","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"The original DFT-raMO was written in MATLAB, but for the sake of performance, reproducibility, and better software architecture, we've rewritten all of the functionality in Julia. You can find the original MATLAB version in our repository.","category":"page"},{"location":"index.html#Publications","page":"Home","title":"Publications","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"Publications outlining the DFT- and Hückel-raMO concepts","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"Yannello, V. J.; Lu, E.; Fredrickson, D. C. At the Limits of Isolobal Bonding: π-Based Covalent Magnetism in Mn2Hg5. Inorg. Chem. 2020, 59 (17), 12304–12313. https://doi.org/10.1021/acs.inorgchem.0c01393.\nYannello, V. J.; Kilduff, B. J.; Fredrickson, D. C. Isolobal Analogies in Intermetallics: The Reversed Approximation MO Approach and Applications to CrGa4- and Ir3Ge7-Type Phases. Inorg. Chem. 2014, 53 (5), 2730–2741. https://doi.org/10.1021/ic4031624.","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"Publications in which raMO is used","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"Kraus, J. D.; Van Buskirk, J. S.; Fredrickson, D. C. The Zintl Concept Applied to Intergrowth Structures: Electron-Hole Matching, Stacking Preferences, and Chemical Pressures in Pd5InAs. Z. Anorg. Allg. Chem. 2023, 649, e202300125. https://doi.org/10.1002/zaac.202300125.\nLim, A.; Fredrickson, D. C. Entropic Control of Bonding, Guided by Chemical Pressure: Phase Transitions and 18-n+m Isomerism of IrIn3. Inorg. Chem. 2023. 62, 27, 10833-10846. https://doi.org/10.1021/acs.inorgchem.3c01496.\nLim, A; Hilleke, K. P.; Fredrickson, D. C. Emergent Transitions: Discord between Electronic and Chemical Pressure Effects in the REAl3 (RE = Sc, Y, Lanthanides) Series. Inorg. Chem. 2023, 62, 4405-4416. https://doi.org/10.1021/acs.inorgchem.2c03393.\nPark, S.-W.; Hosono, H.; Fredrickson, D. C. Cation Clustering in Intermetallics: The Modular Bonding Schemes of CaCu and Ca2Cu. Inorg. Chem. 2019, 58 (15), 10313–10322. https://doi.org/10.1021/acs.inorgchem.9b01486.\nVinokur, A. I.; Fredrickson, D. C. 18-Electron Resonance Structures in the BCC Transition Metals and Their CsCl-Type Derivatives. Inorg. Chem. 2017, 56 (5), 2834–2842. https://doi.org/10.1021/acs.inorgchem.6b02989.\nMiyazaki, K.; Yannello, V. J.; Fredrickson, D. C. Electron-Counting in Intermetallics Made Easy: The 18-n Rule and Isolobal Bonds across the Os–Al System. Zeitschrift für Kristallographie - Crystalline Materials 2017, 232 (7–9), 487–496. https://doi.org/10.1515/zkri-2017-2044.\nEngelkemier, J.; Green, L. M.; McDougald, R. N.; McCandless, G. T.; Chan, J. Y.; Fredrickson, D. C. Putting ScTGa5 (T = Fe, Co, Ni) on the Map: How Electron Counts and Chemical Pressure Shape the Stability Range of the HoCoGa5 Type. Crystal Growth & Design 2016, 16 (9), 5349–5358. https://doi.org/10.1021/acs.cgd.6b00855.\nYannello, V. J.; Fredrickson, D. C. Generality of the 18-n Rule: Intermetallic Structural Chemistry Explained through Isolobal Analogies to Transition Metal Complexes. Inorg. Chem. 2015, 54 (23), 11385–11398. https://doi.org/10.1021/acs.inorgchem.5b02016.\nKilduff, B. J.; Yannello, V. J.; Fredrickson, D. C. Defusing Complexity in Intermetallics: How Covalently Shared Electron Pairs Stabilize the FCC Variant Mo2CuxGa6–x (x ≈ 0.9). Inorg. Chem. 2015, 54 (16), 8103–8110. https://doi.org/10.1021/acs.inorgchem.5b01333. \nYannello, V. J.; Fredrickson, D. C. Orbital Origins of Helices and Magic Electron Counts in the Nowotny Chimney Ladders: The 18 – n Rule and a Path to Incommensurability. Inorg. Chem. 2014, 53 (19), 10627–10631. https://doi.org/10.1021/ic501723n.\nHadler, A. B.; Yannello, V. J.; Bi, W.; Alp, E. E.; Fredrickson, D. C. π-Conjugation in Gd13Fe10C13 and Its Oxycarbide: Unexpected Connections between Complex Carbides and Simple Organic Molecules. J. Am. Chem. Soc. 2014, 136 (34), 12073–12084. https://doi.org/10.1021/ja505868w.","category":"page"}]
}