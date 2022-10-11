# DFTraMO.jl

raMO (reversed approximation Molecular Orbitals) is a Wannier-type analysis that can be used to
reinterpret the results of an electronic structure calculation in terms of localized orbitals,
similar to natural hybrid orbital (NHO) and natural bond orbital (NBO) analysis.

Unlike NBO, raMO follows user-guided targets, which can be centered not just on atoms, but on any
user-selected site. The reconstructions are not dependent on the basis set used to perform the
calculation, which makes this method ideal for solid-state structures.

The original DFT-raMO was written in MATLAB, but for the sake of performance, reproducibility, and
better software architecture, we've rewritten all of the functionality in Julia. You can find the
original MATLAB version in our repository.

# Publications

Have a look at our previous publications on raMO!
1. Yannello, V. J.; Lu, E.; Fredrickson, D. C. At the Limits of Isolobal Bonding: π-Based Covalent Magnetism in Mn2Hg5. _Inorg. Chem._ **2020**, _59_ (17), 12304–12313. https://doi.org/10.1021/acs.inorgchem.0c01393.
2. Park, S.-W.; Hosono, H.; Fredrickson, D. C. Cation Clustering in Intermetallics: The Modular Bonding Schemes of CaCu and Ca2Cu. _Inorg. Chem._ **2019**, _58_ (15), 10313–10322. https://doi.org/10.1021/acs.inorgchem.9b01486.
3. Vinokur, A. I.; Fredrickson, D. C. 18-Electron Resonance Structures in the BCC Transition Metals and Their CsCl-Type Derivatives. _Inorg. Chem._ **2017,** _56_ (5), 2834–2842. https://doi.org/10.1021/acs.inorgchem.6b02989.
4. Miyazaki, K.; Yannello, V. J.; Fredrickson, D. C. Electron-Counting in Intermetallics Made Easy: The 18-n Rule and Isolobal Bonds across the Os–Al System. _Zeitschrift für Kristallographie - Crystalline Materials_ **2017**, _232_ (7–9), 487–496. https://doi.org/10.1515/zkri-2017-2044.
5. Engelkemier, J.; Green, L. M.; McDougald, R. N.; McCandless, G. T.; Chan, J. Y.; Fredrickson, D. C. Putting ScTGa5 (T = Fe, Co, Ni) on the Map: How Electron Counts and Chemical Pressure Shape the Stability Range of the HoCoGa5 Type. _Crystal Growth & Design_ **2016**, _16_ (9), 5349–5358. https://doi.org/10.1021/acs.cgd.6b00855.
6. Yannello, V. J.; Fredrickson, D. C. Generality of the 18-n Rule: Intermetallic Structural Chemistry Explained through Isolobal Analogies to Transition Metal Complexes. _Inorg. Chem._ **2015**, _54_ (23), 11385–11398. https://doi.org/10.1021/acs.inorgchem.5b02016.
7. Kilduff, B. J.; Yannello, V. J.; Fredrickson, D. C. Defusing Complexity in Intermetallics: How Covalently Shared Electron Pairs Stabilize the FCC Variant Mo2CuxGa6–x (x ≈ 0.9). _Inorg. Chem._ **2015**, _54_ (16), 8103–8110. https://doi.org/10.1021/acs.inorgchem.5b01333.
8. Yannello, V. J.; Kilduff, B. J.; Fredrickson, D. C. Isolobal Analogies in Intermetallics: The Reversed Approximation MO Approach and Applications to CrGa4- and Ir3Ge7-Type Phases. _Inorg. Chem._ **2014**, _53_ (5), 2730–2741. https://doi.org/10.1021/ic4031624.
9. Yannello, V. J.; Fredrickson, D. C. Orbital Origins of Helices and Magic Electron Counts in the Nowotny Chimney Ladders: The 18 – n Rule and a Path to Incommensurability. _Inorg. Chem._ **2014,** _53_ (19), 10627–10631. https://doi.org/10.1021/ic501723n.
10. Hadler, A. B.; Yannello, V. J.; Bi, W.; Alp, E. E.; Fredrickson, D. C. π-Conjugation in Gd13Fe10C13 and Its Oxycarbide: Unexpected Connections between Complex Carbides and Simple Organic Molecules. _J. Am. Chem. Soc._ **2014**, _136_ (34), 12073–12084. https://doi.org/10.1021/ja505868w.
