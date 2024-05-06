# DFT calculations

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
