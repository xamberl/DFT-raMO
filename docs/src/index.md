![DFT-raMO.jl logo](assets/dftramologo_500.png)
# DFT-raMO.jl

The DFT-raMO (reversed approximation Molecular Orbital) method is a Wannier-type analysis that can be used to
reinterpret the results of an electronic structure calculation in terms of localized orbitals,
similar to natural hybrid orbital (NHO) and natural bond orbital (NBO) analysis.

Unlike NBO, raMO follows user-guided targets, which can be centered not just on atoms, but on any
user-selected site. The reconstructions are not dependent on the basis set used to perform the
calculation, which makes this method ideal for solid-state structures.

The original DFT-raMO was written in MATLAB, but for the sake of performance, reproducibility, and
better software architecture, we've rewritten all of the functionality in Julia. You can find the
original MATLAB version in our repository.

# Installation
DFT-raMO.jl requires installation of [Julia](https://julialang.org/downloads/) 1.8 or later, available on multiple platforms.
After installation of Julia, the console window (REPL) can be opened via command line by typing the following:
```
julia
```
Next, open the package manager by typing `]`.
This nows allows us to install DFT-raMO.jl with the command:
```
add https://github.com/xamberl/DFT-raMO
```
Wait for the package manager to finish installation, then exit the package manager with the `Backspace` key. If desired, we can exit the Julia REPL by typing:
```
exit()
```

# Publications

Publications outlining the DFT- and Hückel-raMO concepts
* Yannello, V. J.; Lu, E.; Fredrickson, D. C. At the Limits of Isolobal Bonding: π-Based Covalent Magnetism in Mn₂Hg₅. _Inorg. Chem._ **2020**, _59_ (17), 12304–12313. [https://doi.org/10.1021/acs.inorgchem.0c01393](https://doi.org/10.1021/acs.inorgchem.0c01393).
* Yannello, V. J.; Kilduff, B. J.; Fredrickson, D. C. Isolobal Analogies in Intermetallics: The Reversed Approximation MO Approach and Applications to CrGa₄- and Ir₃Ge₇-Type Phases. _Inorg. Chem._ **2014**, _53_ (5), 2730–2741. [https://doi.org/10.1021/ic4031624](https://doi.org/10.1021/ic4031624).

Publications in which raMO is used
* Kraus, J. D.; Van Buskirk, J. S.; Fredrickson, D. C. The Zintl Concept Applied to Intergrowth Structures: Electron-Hole Matching, Stacking Preferences, and Chemical Pressures in Pd₅InAs. _Z. Anorg. Allg. Chem._ **2023**, _649_, e202300125. [https://doi.org/10.1002/zaac.202300125](https://doi.org/10.1002/zaac.202300125).
* Lim, A.; Fredrickson, D. C. Entropic Control of Bonding, Guided by Chemical Pressure: Phase Transitions and 18-*n*+*m* Isomerism of IrIn₃. _Inorg. Chem._ **2023**. _62_, 27, 10833-10846. [https://doi.org/10.1021/acs.inorgchem.3c01496](https://doi.org/10.1021/acs.inorgchem.3c01496).
* Lim, A; Hilleke, K. P.; Fredrickson, D. C. Emergent Transitions: Discord between Electronic and Chemical Pressure Effects in the *RE*Al₃ (*RE* = Sc, Y, Lanthanides) Series. _Inorg. Chem._ **2023**, _62_, 4405-4416. [https://doi.org/10.1021/acs.inorgchem.2c03393](https://doi.org/10.1021/acs.inorgchem.2c03393).
* Park, S.-W.; Hosono, H.; Fredrickson, D. C. Cation Clustering in Intermetallics: The Modular Bonding Schemes of CaCu and Ca₂Cu. _Inorg. Chem._ **2019**, _58_ (15), 10313–10322. [https://doi.org/10.1021/acs.inorgchem.9b01486](https://doi.org/10.1021/acs.inorgchem.9b01486).
* Vinokur, A. I.; Fredrickson, D. C. 18-Electron Resonance Structures in the BCC Transition Metals and Their CsCl-Type Derivatives. _Inorg. Chem._ **2017,** _56_ (5), 2834–2842. [https://doi.org/10.1021/acs.inorgchem.6b02989](https://doi.org/10.1021/acs.inorgchem.6b02989).
* Miyazaki, K.; Yannello, V. J.; Fredrickson, D. C. Electron-Counting in Intermetallics Made Easy: The 18-*n* Rule and Isolobal Bonds across the Os–Al System. _Zeitschrift für Kristallographie - Crystalline Materials_ **2017**, _232_ (7–9), 487–496. [https://doi.org/10.1515/zkri-2017-2044](https://doi.org/10.1515/zkri-2017-2044).
* Engelkemier, J.; Green, L. M.; McDougald, R. N.; McCandless, G. T.; Chan, J. Y.; Fredrickson, D. C. Putting ScTGa₅ (T = Fe, Co, Ni) on the Map: How Electron Counts and Chemical Pressure Shape the Stability Range of the HoCoGa₅ Type. _Crystal Growth & Design_ **2016**, _16_ (9), 5349–5358. [https://doi.org/10.1021/acs.cgd.6b00855](https://doi.org/10.1021/acs.cgd.6b00855).
* Yannello, V. J.; Fredrickson, D. C. Generality of the 18-*n* Rule: Intermetallic Structural Chemistry Explained through Isolobal Analogies to Transition Metal Complexes. _Inorg. Chem._ **2015**, _54_ (23), 11385–11398. [https://doi.org/10.1021/acs.inorgchem.5b02016](https://doi.org/10.1021/acs.inorgchem.5b02016).
* Kilduff, B. J.; Yannello, V. J.; Fredrickson, D. C. Defusing Complexity in Intermetallics: How Covalently Shared Electron Pairs Stabilize the FCC Variant Mo₂CuₓGa₆₋ₓ (x ≈ 0.9). _Inorg. Chem._ **2015**, _54_ (16), 8103–8110. [https://doi.org/10.1021/acs.inorgchem.5b01333](https://doi.org/10.1021/acs.inorgchem.5b01333). 
* Yannello, V. J.; Fredrickson, D. C. Orbital Origins of Helices and Magic Electron Counts in the Nowotny Chimney Ladders: The 18-*n* Rule and a Path to Incommensurability. _Inorg. Chem._ **2014,** _53_ (19), 10627–10631. [https://doi.org/10.1021/ic501723n](https://doi.org/10.1021/ic501723n).
* Hadler, A. B.; Yannello, V. J.; Bi, W.; Alp, E. E.; Fredrickson, D. C. π-Conjugation in Gd₁₃Fe₁₀C₁₃ and Its Oxycarbide: Unexpected Connections between Complex Carbides and Simple Organic Molecules. _J. Am. Chem. Soc._ **2014**, _136_ (34), 12073–12084. [https://doi.org/10.1021/ja505868w](https://doi.org/10.1021/ja505868w).
