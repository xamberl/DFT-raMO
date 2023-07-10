#=
Double semicolons are needed to get the matrix dimensionality correct
Parameters are (in order):
    atomic number
    number of valence electrons
    orbital angular momentum (l)
    principal quantum number (n)
    IP (Energy of the atomic orital)
    Exponents of the STO linear combination
    Scalar coefficients of the STO linear combination
=#
const DFTRAMO_EHT_PARAMS = ehtParams(
    OrbitalParams[
        # Hydrogen
        OrbitalParams(1, 1, 0, 1, -10.0, 6.977549, 0.0, 1.0, 0.0)
        OrbitalParams(1, 1, 1, 2,  -7.0,      6.8, 0.0, 1.0, 0.0)
        OrbitalParams(1, 1, 2, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(1, 1, 3, 0,   0.0,      0.0, 0.0, 0.0, 0.0);;
        # Helium
        OrbitalParams(2, 2, 0, 1, -10.0, 1.620030, 0.0, 1.0, 0.0)
        OrbitalParams(2, 2, 1, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(2, 2, 2, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(2, 2, 3, 0,   0.0,      0.0, 0.0, 0.0, 0.0);;
        # Lithium
        OrbitalParams(3, 1, 0, 2, -10.0, 0.510296, 0.0, 1.0, 0.0)
        OrbitalParams(3, 1, 1, 2,  -7.0, 0.565082, 0.0, 1.0, 0.0)
        OrbitalParams(3, 1, 2, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(3, 1, 3, 0,   0.0,      0.0, 0.0, 0.0, 0.0);;
        # Beryllium
        OrbitalParams(4, 2, 0, 2, -10.0, 0.954418, 0.0, 1.0, 0.0)
        OrbitalParams(4, 2, 1, 2,  -7.0, 0.922233, 0.0, 1.0, 0.0)
        OrbitalParams(4, 2, 2, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(4, 2, 3, 0,   0.0,      0.0, 0.0, 0.0, 0.0);;
        # Boron
        OrbitalParams(5, 3, 0, 2, -10.0, 1.331883, 0.0, 1.0, 0.0)
        OrbitalParams(5, 3, 1, 2,  -7.0, 1.270451, 0.0, 1.0, 0.0)
        OrbitalParams(5, 3, 2, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(5, 3, 3, 0,   0.0,      0.0, 0.0, 0.0, 0.0);;
        # Carbon
        OrbitalParams(6, 4, 0, 2, -10.0, 1.676061, 0.0, 1.0, 0.0)
        OrbitalParams(6, 4, 1, 2,  -7.0, 1.605091, 0.0, 1.0, 0.0)
        OrbitalParams(6, 4, 2, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(6, 4, 3, 0,   0.0,      0.0, 0.0, 0.0, 0.0);;
        # Nitrogen
        OrbitalParams(7, 5, 0, 2, -10.0, 2.008054, 0.0, 1.0, 0.0)
        OrbitalParams(7, 5, 1, 2,  -7.0, 1.825051, 0.0, 1.0, 0.0)
        OrbitalParams(7, 5, 2, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(7, 5, 3, 0,   0.0,      0.0, 0.0, 0.0, 0.0);;
        # Oxygen
        OrbitalParams(8, 6, 0, 2, -10.0, 2.281453, 0.0, 1.0, 0.0)
        OrbitalParams(8, 6, 1, 2,  -7.0, 1.991072, 0.0, 1.0, 0.0)
        OrbitalParams(8, 6, 2, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(8, 6, 3, 0,   0.0,      0.0, 0.0, 0.0, 0.0);;
        # Fluorine
        OrbitalParams(9, 7, 0, 2, -10.0, 2.530891, 0.0, 1.0, 0.0)
        OrbitalParams(9, 7, 1, 2,  -7.0, 2.147303, 0.0, 1.0, 0.0)
        OrbitalParams(9, 7, 2, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(9, 7, 3, 0,   0.0,      0.0, 0.0, 0.0, 0.0);;
        # Neon
        OrbitalParams(10, 8, 0, 2, -10.0, 2.656822, 0.0, 1.0, 0.0)
        OrbitalParams(10, 8, 1, 2,  -7.0, 2.144320, 0.0, 1.0, 0.0)
        OrbitalParams(10, 8, 2, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(10, 8, 3, 0,   0.0,      0.0, 0.0, 0.0, 0.0);;
        # Sodium
        OrbitalParams(11, 1, 0, 3, -10.0, 0.943660, 0.0, 1.0, 0.0)
        OrbitalParams(11, 1, 1, 3,  -7.0, 0.591956, 0.0, 1.0, 0.0)
        OrbitalParams(11, 1, 2, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(11, 1, 3, 0,   0.0,      0.0, 0.0, 0.0, 0.0);;
        # Magnesium
        OrbitalParams(12, 2, 0, 3, -10.0, 1.120316, 0.0, 1.0, 0.0)
        OrbitalParams(12, 2, 1, 3,  -7.0, 0.799968, 0.0, 1.0, 0.0)
        OrbitalParams(12, 2, 2, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(12, 2, 3, 0,   0.0,      0.0, 0.0, 0.0, 0.0);;
        # Aluminum
        OrbitalParams(13, 3, 0, 3, -10.0, 1.410446, 0.0, 1.0, 0.0)
        OrbitalParams(13, 3, 1, 3,  -7.0, 1.119116, 0.0, 1.0, 0.0)
        OrbitalParams(13, 3, 2, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(13, 3, 3, 0,   0.0,      0.0, 0.0, 0.0, 0.0);;
        # Silicon
        OrbitalParams(14, 4, 0, 3, -10.0, 1.703463, 0.0, 1.0, 0.0)
        OrbitalParams(14, 4, 1, 3,  -7.0, 1.410385, 0.0, 1.0, 0.0)
        OrbitalParams(14, 4, 2, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(14, 4, 3, 0,   0.0,      0.0, 0.0, 0.0, 0.0);;
        # Phosphorus
        OrbitalParams(15, 5, 0, 3, -10.0, 1.969166, 0.0, 1.0, 0.0)
        OrbitalParams(15, 5, 1, 3,  -7.0, 1.662988, 0.0, 1.0, 0.0)
        OrbitalParams(15, 5, 2, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(15, 5, 3, 0,   0.0,      0.0, 0.0, 0.0, 0.0);;
        # Sulfur
        OrbitalParams(16, 6, 0, 3, -10.0, 2.197697, 0.0, 1.0, 0.0)
        OrbitalParams(16, 6, 1, 3,  -7.0, 1.871105, 0.0, 1.0, 0.0)
        OrbitalParams(16, 6, 2, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(16, 6, 3, 0,   0.0,      0.0, 0.0, 0.0, 0.0);;
        # Chlorine
        OrbitalParams(17, 7, 0, 3, -10.0, 2.381834, 0.0, 1.0, 0.0)
        OrbitalParams(17, 7, 1, 3,  -7.0, 2.031865, 0.0, 1.0, 0.0)
        OrbitalParams(17, 7, 2, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(17, 7, 3, 0,   0.0,      0.0, 0.0, 0.0, 0.0);;
        # Argon
        OrbitalParams(18, 8, 0, 3, -10.0, 2.555607, 0.0, 1.0, 0.0)
        OrbitalParams(18, 8, 1, 3,  -7.0, 2.203441, 0.0, 1.0, 0.0)
        OrbitalParams(18, 8, 2, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(18, 8, 3, 0,   0.0,      0.0, 0.0, 0.0, 0.0);;
        # Potassium
        OrbitalParams(19, 1, 0, 4, -10.0, 0.975319,      0.0,        1.0,       0.0)
        OrbitalParams(19, 1, 1, 4,  -3.0, 0.221807,      0.0,        1.0,       0.0)
        OrbitalParams(19, 1, 2, 3,  -7.0, 0.464396, 1.809022, 184.919580, 20.631909) # sus
        OrbitalParams(19, 1, 3, 0,   0.0,      0.0,      0.0,        0.0,       0.0);;
        # Calcium
        OrbitalParams(20, 2, 0, 4, -10.0, 1.181317,      0.0,      1.0,      0.0)
        OrbitalParams(20, 2, 1, 4,  -3.0, 0.941368,      0.0,      1.0,      0.0)
        OrbitalParams(20, 2, 2, 3,  -7.0, 0.846255, 1.315737, 0.439826, 0.956826)
        OrbitalParams(20, 2, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Scandium
        OrbitalParams(21, 3, 0, 4,  -7.0, 1.263174,      0.0,      1.0,      0.0)
        OrbitalParams(21, 3, 1, 4,  -3.0, 0.741697,      0.0,      1.0,      0.0)
        OrbitalParams(21, 3, 2, 3, -10.0, 0.440416, 1.518473, 1.272698, 1.947687)
        OrbitalParams(21, 3, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Titanium
        OrbitalParams(22, 4, 0, 4,  -7.0, 1.342962,      0.0,      1.0,     0.0)
        OrbitalParams(22, 4, 1, 4,  -3.0, 0.821780,      0.0,      1.0,     0.0)
        OrbitalParams(22, 4, 2, 3, -10.0, 0.300286, 1.675023, 1.962965, 1.32482)
        OrbitalParams(22, 4, 3, 0,   0.0,      0.0,      0.0,      0.0,     0.0);;
        # Vanadium
        OrbitalParams(23, 5, 0, 4,  -7.0, 1.425320,      0.0,      1.0,      0.0)
        OrbitalParams(23, 5, 1, 4,  -3.0, 1.014350,      0.0,      1.0,      0.0)
        OrbitalParams(23, 5, 2, 3, -10.0, 0.277856, 1.749907, 1.263332, 0.931017)
        OrbitalParams(23, 5, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Chromium
        OrbitalParams(24, 6, 0, 4,  -7.0, 1.492470,      0.0,      1.0,      0.0)
        OrbitalParams(24, 6, 1, 4,  -3.0, 0.961316,      0.0,      1.0,      0.0)
        OrbitalParams(24, 6, 2, 3, -10.0, 0.295693, 1.952856, 2.934024, 1.707039)
        OrbitalParams(24, 6, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Manganese
        OrbitalParams(25, 7, 0, 4,  -7.0, 1.567448,      0.0,      1.0,      0.0)
        OrbitalParams(25, 7, 1, 4,  -3.0, 1.118844,      0.0,      1.0,      0.0)
        OrbitalParams(25, 7, 2, 3, -10.0, 0.270335, 2.091848, 2.971582, 1.590746)
        OrbitalParams(25, 7, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Iron
        OrbitalParams(26, 8, 0, 4,  -7.0, 1.628480,      0.0,      1.0,      0.0)
        OrbitalParams(26, 8, 1, 4,  -3.0, 1.185066,      0.0,      1.0,      0.0)
        OrbitalParams(26, 8, 2, 3, -10.0, 0.268476, 2.170869, 3.352271, 1.338474)
        OrbitalParams(26, 8, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Cobalt
        OrbitalParams(27, 9, 0, 4,  -7.0, 1.698666,      0.0,      1.0,      0.0)
        OrbitalParams(27, 9, 1, 4,  -3.0, 1.243328,      0.0,      1.0,      0.0)
        OrbitalParams(27, 9, 2, 3, -10.0, 0.284904, 2.162554, 2.547326, 1.747923)
        OrbitalParams(27, 9, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Nickel
        OrbitalParams(28, 10, 0, 4,  -7.0, 1.753755,      0.0,      1.0,     0.0)
        OrbitalParams(28, 10, 1, 4,  -3.0, 1.286072,      0.0,      1.0,     0.0)
        OrbitalParams(28, 10, 2, 3, -10.0, 0.229392, 2.186438, 3.036308, 1.37272)
        OrbitalParams(28, 10, 3, 0,   0.0,      0.0,      0.0,      0.0,     0.0);;
        # Copper
        OrbitalParams(29, 1, 0, 4,  -7.0, 1.828572,      0.0,      1.0,      0.0)
        OrbitalParams(29, 1, 1, 4,  -3.0, 1.345442,      0.0,      1.0,      0.0)
        OrbitalParams(29, 1, 2, 3, -10.0, 0.088592, 2.150036, 4.079091, 1.090018)
        OrbitalParams(29, 1, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Zinc
        OrbitalParams(30, 2, 0, 4,  -7.0, 1.938283,      0.0,      1.0,      0.0)
        OrbitalParams(30, 2, 1, 4,  -3.0, 1.466064,      0.0,      1.0,      0.0)
        OrbitalParams(30, 2, 2, 3, -10.0, 0.527206, 2.315908, 0.000008, 2.196275)
        OrbitalParams(30, 2, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Gallium
        OrbitalParams(31, 3, 0, 4, -10.0, 2.011001, 0.0, 1.0, 0.0)
        OrbitalParams(31, 3, 1, 4,  -7.0, 1.617839, 0.0, 1.0, 0.0)
        OrbitalParams(31, 3, 2, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(31, 3, 3, 0,   0.0,      0.0, 0.0, 0.0, 0.0);;
        # Germanium
        OrbitalParams(32, 4, 0, 4, -10.0, 2.263211, 0.0, 1.0, 0.0)
        OrbitalParams(32, 4, 1, 4,  -7.0, 1.859698, 0.0, 1.0, 0.0)
        OrbitalParams(32, 4, 2, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(32, 4, 3, 0,   0.0,      0.0, 0.0, 0.0, 0.0);;
        # Arsenic
        OrbitalParams(33, 5, 0, 4, -10.0, 2.462542, 0.0, 1.0, 0.0)
        OrbitalParams(33, 5, 1, 4,  -7.0, 2.052448, 0.0, 1.0, 0.0)
        OrbitalParams(33, 5, 2, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(33, 5, 3, 0,   0.0,      0.0, 0.0, 0.0, 0.0);;
        # Selenium
        OrbitalParams(34, 6, 0, 4, -10.0, 2.603851, 0.0, 1.0, 0.0)
        OrbitalParams(34, 6, 1, 4,  -7.0, 2.187979, 0.0, 1.0, 0.0)
        OrbitalParams(34, 6, 2, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(34, 6, 3, 0,   0.0,      0.0, 0.0, 0.0, 0.0);;
        # Bromine
        OrbitalParams(35, 7, 0, 4, -10.0, 2.747563, 0.0, 1.0, 0.0)
        OrbitalParams(35, 7, 1, 4,  -7.0, 2.321994, 0.0, 1.0, 0.0)
        OrbitalParams(35, 7, 2, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(35, 7, 3, 0,   0.0,      0.0, 0.0, 0.0, 0.0);;
        # Krypton
        OrbitalParams(36, 8, 0, 4, -10.0, 2.801036, 0.0, 1.0, 0.0)
        OrbitalParams(36, 8, 1, 4,  -7.0, 2.375968, 0.0, 1.0, 0.0)
        OrbitalParams(36, 8, 2, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(36, 8, 3, 0,   0.0,      0.0, 0.0, 0.0, 0.0);;
        # Rubidium
        OrbitalParams(37, 1, 0, 5, -10.0, 1.220789,      0.0,      1.0,      0.0)
        OrbitalParams(37, 1, 1, 5,  -3.0, 0.374380,      0.0,      1.0,      0.0)
        OrbitalParams(37, 1, 2, 4,  -7.0, 0.496336, 1.710271, 3.996439, 0.822399)
        OrbitalParams(37, 1, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Strontium
        OrbitalParams(38, 2, 0, 5, -10.0, 1.340750,      0.0,      1.0,      0.0)
        OrbitalParams(38, 2, 1, 5,  -3.0, 0.517207,      0.0,      1.0,      0.0)
        OrbitalParams(38, 2, 2, 4,  -7.0, 0.635090, 2.213849, 2.476703, 0.969176)
        OrbitalParams(38, 2, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Yttrium
        OrbitalParams(39, 3, 0, 5,  -7.0, 1.460711,      0.0,      1.0,      0.0)
        OrbitalParams(39, 3, 1, 5,  -3.0, 0.478628,      0.0,      1.0,      0.0)
        OrbitalParams(39, 3, 2, 4, -10.0, 0.555035, 2.183579, 2.107518, 0.855612)
        OrbitalParams(39, 3, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Zirconium
        OrbitalParams(40, 4, 0, 5,  -7.0, 1.580671,      0.0,     1.0,      0.0)
        OrbitalParams(40, 4, 1, 5,  -3.0, 0.988008,      0.0,     1.0,      0.0)
        OrbitalParams(40, 4, 2, 4, -10.0, 0.527534, 2.003386, 1.97010, 1.498781)
        OrbitalParams(40, 4, 3, 0,   0.0,      0.0,      0.0,     0.0,      0.0);;
        # Niobium
        OrbitalParams(41, 5, 0, 5,  -7.0, 1.670563,      0.0,      1.0,      0.0)
        OrbitalParams(41, 5, 1, 5,  -3.0, 1.035008,      0.0,      1.0,      0.0)
        OrbitalParams(41, 5, 2, 4, -10.0, 0.574200, 2.315978, 1.790895, 1.210563)
        OrbitalParams(41, 5, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Molybdenum
        OrbitalParams(42, 6, 0, 5,  -7.0, 1.733234,      0.0,      1.0,      0.0)
        OrbitalParams(42, 6, 1, 5,  -3.0, 1.625703,      0.0,      1.0,      0.0)
        OrbitalParams(42, 6, 2, 4, -10.0, 0.536506, 2.422384, 2.035252, 1.237972)
        OrbitalParams(42, 6, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Technetium
        OrbitalParams(43, 7, 0, 5,  -7.0, 1.792823,      0.0,      1.0,     0.0)
        OrbitalParams(43, 7, 1, 5,  -3.0, 1.718796,      0.0,      1.0,     0.0)
        OrbitalParams(43, 7, 2, 4, -10.0, 0.565373, 2.488641, 2.081664, 1.65479)
        OrbitalParams(43, 7, 3, 0,   0.0,      0.0,      0.0,      0.0,     0.0);;
        # Ruthenium
        OrbitalParams(44, 8, 0, 5,  -7.0, 1.875697,      0.0,      1.0,     0.0)
        OrbitalParams(44, 8, 1, 5,  -3.0, 1.387973,      0.0,      1.0,     0.0)
        OrbitalParams(44, 8, 2, 4, -10.0, 0.786233, 2.520292, 1.273693, 2.55724)
        OrbitalParams(44, 8, 3, 0,   0.0,      0.0,      0.0,      0.0,     0.0);;
        # Rhodium
        OrbitalParams(45, 9, 0, 5,  -7.0, 1.940834,      0.0,      1.0,      0.0)
        OrbitalParams(45, 9, 1, 5,  -3.0, 1.448122,      0.0,      1.0,      0.0)
        OrbitalParams(45, 9, 2, 4, -10.0, 0.595925, 2.581448, 1.510618, 1.899981)
        OrbitalParams(45, 9, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Palladium
        OrbitalParams(46, 10, 0, 5,  -7.0, 2.025704,     0.0,      1.0,      0.0)
        OrbitalParams(46, 10, 1, 5,  -3.0, 1.500623,     0.0,      1.0,      0.0)
        OrbitalParams(46, 10, 2, 4, -10.0, 0.670094, 2.70693, 0.877034, 1.482471)
        OrbitalParams(46, 10, 3, 0,   0.0,      0.0,     0.0,      0.0,      0.0);;
        # Silver
        OrbitalParams(47, 11, 0, 5,  -7.0, 2.079598,      0.0,      1.0,      0.0)
        OrbitalParams(47, 11, 1, 5,  -3.0, 1.559465,      0.0,      1.0,      0.0)
        OrbitalParams(47, 11, 2, 4, -10.0, 1.252589, 2.788983, 0.317566, 1.540176)
        OrbitalParams(47, 11, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Cadmium
        OrbitalParams(48, 12, 0, 5,  -7.0, 2.165696,      0.0,     1.0,      0.0)
        OrbitalParams(48, 12, 1, 5,  -3.0, 1.675763,      0.0,     1.0,      0.0)
        OrbitalParams(48, 12, 2, 4, -10.0, 0.422275, 2.952111, 2.20011, 1.048571)
        OrbitalParams(48, 12, 3, 0,   0.0,      0.0,      0.0,     0.0,      0.0);;
        # Indium
        OrbitalParams(49, 3, 0, 5, -10.0, 2.118550, 0.0, 1.0, 0.0)
        OrbitalParams(49, 3, 1, 5,  -7.0, 1.740947, 0.0, 1.0, 0.0)
        OrbitalParams(49, 3, 2, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(49, 3, 3, 0,   0.0,      0.0, 0.0, 0.0, 0.0);;
        # Tin
        OrbitalParams(50, 4, 0, 5, -10.0, 2.293511, 0.0, 1.0, 0.0)
        OrbitalParams(50, 4, 1, 5,  -7.0, 1.909329, 0.0, 1.0, 0.0)
        OrbitalParams(50, 4, 2, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(50, 4, 3, 0,   0.0,      0.0, 0.0, 0.0, 0.0);;
        # Antimony
        OrbitalParams(51, 5, 0, 5, -10.0, 2.606718, 0.0, 1.0, 0.0)
        OrbitalParams(51, 5, 1, 5,  -7.0, 2.189100, 0.0, 1.0, 0.0)
        OrbitalParams(51, 5, 2, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(51, 5, 3, 0,   0.0,      0.0, 0.0, 0.0, 0.0);;
        # Tellurium
        OrbitalParams(52, 6, 0, 5, -10.0, 2.723071, 0.0, 1.0, 0.0)
        OrbitalParams(52, 6, 1, 5,  -7.0, 2.308011, 0.0, 1.0, 0.0)
        OrbitalParams(52, 6, 2, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(52, 6, 3, 0,   0.0,      0.0, 0.0, 0.0, 0.0);;
        # Iodine
        OrbitalParams(53, 7, 0, 5, -10.0, 2.854851, 0.0, 1.0, 0.0)
        OrbitalParams(53, 7, 1, 5,  -7.0, 2.435386, 0.0, 1.0, 0.0)
        OrbitalParams(53, 7, 2, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(53, 7, 3, 0,   0.0,      0.0, 0.0, 0.0, 0.0);;
        # Xenon
        OrbitalParams(54, 8, 0, 5, -10.0, 2.888075, 0.0, 1.0, 0.0)
        OrbitalParams(54, 8, 1, 5,  -7.0, 2.488795, 0.0, 1.0, 0.0)
        OrbitalParams(54, 0, 2, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(54, 0, 3, 0,   0.0,      0.0, 0.0, 0.0, 0.0);;
        # Cesium
        OrbitalParams(55, 1, 0, 6, -10.0, 1.582495, 0.0, 1.0, 0.0)
        OrbitalParams(55, 1, 1, 6,  -3.0, 2.088303, 0.0, 1.0, 0.0)
        OrbitalParams(55, 0, 2, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(55, 0, 3, 0,   0.0,      0.0, 0.0, 0.0, 0.0);;
        # Barium
        OrbitalParams(56, 2, 0, 6, -10.0, 1.583879,      0.0,      1.0,      0.0)
        OrbitalParams(56, 2, 1, 6,  -3.0, 2.171974,      0.0,      1.0,      0.0)
        OrbitalParams(56, 2, 2, 5,  -7.0, 0.769348, 2.331612, 1.606572, 0.737765)
        OrbitalParams(56, 0, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Lanthanum
        OrbitalParams(57, 3, 0, 6,  -7.0, 1.585264,      0.0,      1.0,      0.0)
        OrbitalParams(57, 3, 1, 6,  -3.0, 2.255645,      0.0,      1.0,      0.0)
        OrbitalParams(57, 3, 2, 5, -10.0, 0.793984, 2.499550, 1.459186, 0.747428)
        OrbitalParams(57, 3, 3, 4,  -5.0, 0.204058, 2.315306, 2.099396, 1.513016);;
        # Cerium
        OrbitalParams(58, 4, 0, 6,  -5.0, 1.586648,      0.0,      1.0,      0.0)
        OrbitalParams(58, 4, 1, 6,  -3.0, 1.008433,      0.0,      1.0,      0.0)
        OrbitalParams(58, 4, 2, 5, -10.0, 0.861776, 2.185814, 1.637674, 1.517998)
        OrbitalParams(58, 4, 3, 4,  -7.0, 0.234407, 1.866260, 4.196012, 1.181565);;
        # Praseodymium
        OrbitalParams(59, 5, 0, 6,  -7.0, 1.611790,      0.0,      1.0,      0.0)
        OrbitalParams(59, 5, 1, 6,  -3.0, 1.562551,      0.0,      1.0,      0.0)
        OrbitalParams(59, 5, 2, 5, -10.0, 0.745038, 2.504876, 2.200758, 1.055402)
        OrbitalParams(59, 5, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Neodymium
        OrbitalParams(60, 6, 0, 6,  -7.0, 1.636933,      0.0,      1.0,      0.0)
        OrbitalParams(60, 6, 1, 6,  -3.0, 1.618752,      0.0,      1.0,      0.0)
        OrbitalParams(60, 6, 2, 5, -10.0, 0.784601, 2.524242, 1.725767, 0.951561)
        OrbitalParams(60, 6, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Promethium
        OrbitalParams(61, 7, 0, 6,  -7.0, 1.662075,      0.0,      1.0,     0.0)
        OrbitalParams(61, 7, 1, 6,  -3.0, 0.782919,      0.0,      1.0,     0.0)
        OrbitalParams(61, 7, 2, 5, -10.0, 0.761147, 2.528618, 1.574568, 0.80526)
        OrbitalParams(61, 7, 3, 0,   0.0,      0.0,      0.0,      0.0,     0.0);;
        # Samarium
        OrbitalParams(62, 8, 0, 6,  -7.0, 1.687217,      0.0,      1.0,      0.0)
        OrbitalParams(62, 8, 1, 6,  -3.0, 1.625054,      0.0,      1.0,      0.0)
        OrbitalParams(62, 8, 2, 5, -10.0, 0.829042, 2.492654, 2.149629, 0.997777)
        OrbitalParams(62, 8, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Europium
        OrbitalParams(63, 9, 0, 6,  -7.0, 1.712360,      0.0,      1.0,      0.0)
        OrbitalParams(63, 9, 1, 6,  -3.0, 0.513608,      0.0,      1.0,      0.0)
        OrbitalParams(63, 9, 2, 5, -10.0, 0.681006, 2.157792, 1.767884, 1.196325)
        OrbitalParams(63, 9, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Gadolinium
        OrbitalParams(64, 10, 0, 6,  -7.0, 1.829944,      0.0,      1.0,      0.0)
        OrbitalParams(64, 10, 1, 6,  -3.0, 0.627643,      0.0,      1.0,      0.0)
        OrbitalParams(64, 10, 2, 5, -10.0, 0.593889, 2.372038, 1.997234, 0.667022)
        OrbitalParams(64, 10, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Terbium
        OrbitalParams(65, 11, 0, 6,  -7.0, 1.856014,      0.0,      1.0,      0.0)
        OrbitalParams(65, 11, 1, 6,  -3.0, 0.620841,      0.0,      1.0,      0.0)
        OrbitalParams(65, 11, 2, 5, -10.0, 0.886617, 2.454925, 1.495028, 1.297177)
        OrbitalParams(65, 11, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Dysprosium
        OrbitalParams(66, 12, 0, 6,  -7.0, 1.879440,      0.0,      1.0,      0.0)
        OrbitalParams(66, 12, 1, 6,  -3.0, 0.606429,      0.0,      1.0,      0.0)
        OrbitalParams(66, 12, 2, 5, -10.0, 1.003398, 2.497901, 1.072132, 1.087827)
        OrbitalParams(66, 12, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Holmium
        OrbitalParams(67, 13, 0, 6,  -7.0, 1.903055,      0.0,      1.0,      0.0)
        OrbitalParams(67, 13, 1, 6,  -3.0, 0.589964,      0.0,      1.0,      0.0)
        OrbitalParams(67, 13, 2, 5, -10.0, 1.196878, 2.593773, 0.834799, 0.905168)
        OrbitalParams(67, 13, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Erbium
        OrbitalParams(68, 14, 0, 6, -7.0, 1.923472,      0.0,      1.0,      0.0)
        OrbitalParams(68, 14, 1, 6, -3.0, 0.578557,      0.0,      1.0,      0.0)
        OrbitalParams(68, 14, 2, 5, -10.0, 0.57334, 2.353352, 0.662178, 0.195576)
        OrbitalParams(68, 14, 3, 0,   0.0,     0.0,      0.0,      0.0,      0.0);;
        # Thulium
        OrbitalParams(69, 15, 0, 6,  -7.0, 1.944693,      0.0,      1.0,      0.0)
        OrbitalParams(69, 15, 1, 6,  -3.0, 0.565625,      0.0,      1.0,      0.0)
        OrbitalParams(69, 15, 2, 5, -10.0, 0.577834, 2.345774, 1.067049, 0.322703)
        OrbitalParams(69, 15, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Ytterbium
        OrbitalParams(70, 16, 0, 6,  -7.0, 1.898618,      0.0,      1.0,      0.0)
        OrbitalParams(70, 16, 1, 6,  -3.0, 0.566010,      0.0,      1.0,      0.0)
        OrbitalParams(70, 16, 2, 5, -10.0, 0.911065, 2.467122, 1.297333, 0.469341)
        OrbitalParams(70, 16, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Lutetium
        OrbitalParams(71, 17, 0, 6,  -7.0, 2.026001,      0.0,      1.0,      0.0)
        OrbitalParams(71, 17, 1, 6,  -3.0, 0.621305,      0.0,      1.0,      0.0)
        OrbitalParams(71, 17, 2, 5, -10.0, 0.943219, 2.477175, 1.367228, 1.177759)
        OrbitalParams(71, 17, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Hafnium
        OrbitalParams(72, 4, 0, 6,  -7.0, 2.058378,      0.0,      1.0,      0.0)
        OrbitalParams(72, 4, 1, 6,  -3.0, 1.407755,      0.0,      1.0,      0.0)
        OrbitalParams(72, 4, 2, 5, -10.0, 1.315594, 2.887927, 0.363474, 0.377102)
        OrbitalParams(72, 4, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Tantalum
        OrbitalParams(73, 5, 0, 6,  -7.0, 2.120564,      0.0,      1.0,      0.0)
        OrbitalParams(73, 5, 1, 6,  -3.0, 1.654846,      0.0,      1.0,      0.0)
        OrbitalParams(73, 5, 2, 5, -10.0, 0.818876, 2.655920, 2.034779, 1.557278)
        OrbitalParams(73, 5, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Tungsten
        OrbitalParams(74, 6, 0, 6,  -7.0, 2.220140,      0.0,      1.0,      0.0)
        OrbitalParams(74, 6, 1, 6,  -3.0, 1.910224,      0.0,      1.0,      0.0)
        OrbitalParams(74, 6, 2, 5, -10.0, 0.753708, 2.794688, 1.716222, 1.020458)
        OrbitalParams(74, 6, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Rhenium
        OrbitalParams(75, 7, 0, 6,  -7.0, 2.321055,      0.0,      1.0,      0.0)
        OrbitalParams(75, 7, 1, 6,  -3.0, 1.969623,      0.0,      1.0,      0.0)
        OrbitalParams(75, 7, 2, 5, -10.0, 1.111103, 3.022930, 0.772989, 1.034148)
        OrbitalParams(75, 7, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Osmium
        OrbitalParams(76, 8, 0, 6,  -7.0, 2.372426,      0.0,      1.0,      0.0)
        OrbitalParams(76, 8, 1, 6,  -3.0, 1.793455,      0.0,      1.0,      0.0)
        OrbitalParams(76, 8, 2, 5, -10.0, 0.693757, 2.969067, 2.071311, 1.026248)
        OrbitalParams(76, 8, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Iridium
        OrbitalParams(77, 9, 0, 6,  -7.0, 2.427836,      0.0,      1.0,      0.0)
        OrbitalParams(77, 9, 1, 6,  -3.0, 1.826549,      0.0,      1.0,      0.0)
        OrbitalParams(77, 9, 2, 5, -10.0, 0.643802, 2.945609, 2.232089, 1.076656)
        OrbitalParams(77, 9, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Platinum
        OrbitalParams(78, 10, 0, 6,  -7.0, 2.527887,      0.0,      1.0,      0.0)
        OrbitalParams(78, 10, 1, 6,  -3.0, 1.894722,      0.0,      1.0,      0.0)
        OrbitalParams(78, 10, 2, 5, -10.0, 0.714973, 3.099921, 1.576716, 0.956136)
        OrbitalParams(78, 10, 3, 0,   0.0,      0.0,      0.0,      0.0,      0.0);;
        # Gold
        OrbitalParams(79, 11, 0, 6,  -7.0, 2.571008,      0.0,     1.0,      0.0)
        OrbitalParams(79, 11, 1, 6,  -3.0, 1.935752,      0.0,     1.0,      0.0)
        OrbitalParams(79, 11, 2, 5, -10.0, 0.741844, 3.135661, 1.16673, 0.889816)
        OrbitalParams(79, 11, 3, 0,   0.0,      0.0,      0.0,     0.0,      0.0);; 
        # Mercury
        OrbitalParams(80, 12, 0, 6,  -7.0, 2.650556,      0.0,     1.0,      0.0)
        OrbitalParams(80, 12, 1, 6,  -3.0, 2.028762,      0.0,     1.0,      0.0)
        OrbitalParams(80, 12, 2, 5, -10.0, 0.624733, 3.208542, 1.78873, 0.975926)
        OrbitalParams(80, 12, 3, 0,   0.0,      0.0,      0.0,     0.0,      0.0);;
        # Thallium
        OrbitalParams(81, 3, 0, 6, -10.0, 2.470503, 0.0, 1.0, 0.0)
        OrbitalParams(81, 3, 1, 6,  -7.0, 1.991451, 0.0, 1.0, 0.0)
        OrbitalParams(81, 3, 2, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(81, 3, 3, 0,   0.0,      0.0, 0.0, 0.0, 0.0);;
        # Lead
        OrbitalParams(82, 4, 0, 6, -10.0, 2.656227, 0.0, 1.0, 0.0)
        OrbitalParams(82, 4, 1, 6,  -7.0, 2.164111, 0.0, 1.0, 0.0)
        OrbitalParams(82, 4, 2, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(82, 4, 3, 0,   0.0,      0.0, 0.0, 0.0, 0.0);;
        # Bismuth
        OrbitalParams(83, 5, 0, 6, -10.0, 2.847450, 0.0, 1.0, 0.0)
        OrbitalParams(83, 5, 1, 6,  -7.0, 2.339514, 0.0, 1.0, 0.0)
        OrbitalParams(83, 5, 2, 0,   0.0,      0.0, 0.0, 0.0, 0.0)
        OrbitalParams(83, 5, 3, 0,   0.0,      0.0, 0.0, 0.0, 0.0);;
        # Polonium
        OrbitalParams(84, 6, 0, 6, 0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(84, 6, 1, 6, 0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(84, 6, 2, 0, 0.0, 0.0, 0.0, 0.0, 0.0)
        OrbitalParams(84, 6, 3, 0, 0.0, 0.0, 0.0, 0.0, 0.0);;
        # Astatine
        OrbitalParams(85, 7, 0, 6, 0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(85, 7, 1, 6, 0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(85, 7, 2, 0, 0.0, 0.0, 0.0, 0.0, 0.0)
        OrbitalParams(85, 7, 3, 0, 0.0, 0.0, 0.0, 0.0, 0.0);;
        # Radon
        OrbitalParams(86, 8, 0, 6, 0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(86, 8, 1, 6, 0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(86, 8, 2, 0, 0.0, 0.0, 0.0, 0.0, 0.0)
        OrbitalParams(86, 8, 3, 0, 0.0, 0.0, 0.0, 0.0, 0.0);;
        # Francium
        OrbitalParams(87, 1, 0, 7, 0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(87, 1, 1, 7, 0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(87, 1, 2, 0, 0.0, 0.0, 0.0, 0.0, 0.0)
        OrbitalParams(87, 1, 3, 0, 0.0, 0.0, 0.0, 0.0, 0.0);;
        # Radium
        OrbitalParams(88, 2, 0, 7, 0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(88, 2, 1, 7, 0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(88, 2, 2, 0, 0.0, 0.0, 0.0, 0.0, 0.0)
        OrbitalParams(88, 2, 3, 0, 0.0, 0.0, 0.0, 0.0, 0.0);;
        # Actinium
        OrbitalParams(89, 3, 0, 7,  -5.0, 1.920034,      0.0,      1.0,      0.0)
        OrbitalParams(89, 3, 1, 7,  -3.0, 0.902683,      0.0,      1.0,      0.0)
        OrbitalParams(89, 3, 2, 6, -10.0, 0.928908, 2.439746, 1.576864, 1.075456)
        OrbitalParams(89, 3, 3, 5,  -7.0, 0.559570, 2.232725, 2.167963, 1.246645);;
        # Thorium
        OrbitalParams(90, 0, 0, 7,  -5.0, 1.972647,      0.0,      1.0,      0.0)
        OrbitalParams(90, 0, 1, 7,  -3.0, 1.729200,      0.0,      1.0,      0.0)
        OrbitalParams(90, 0, 2, 6, -10.0, 0.864878, 2.440685, 2.156274, 1.349211)
        OrbitalParams(90, 0, 3, 5,  -7.0, 0.459461, 2.291643, 3.145788, 1.331615);;
        # Protactinium
        OrbitalParams(91, 0, 0, 7, -5.0, 2.084153, 0.0, 1.0, 0.0)
        OrbitalParams(91, 0, 1, 7, -3.0, 1.742004, 0.0, 1.0, 0.0)
        OrbitalParams(91, 0, 2, 6, -10.0, 0.783779, 2.582146, 2.381364, 1.054682)
        OrbitalParams(91, 0, 3, 5, -7.0, 0.382723, 2.529159, 3.648382, 1.301798);;
        # Uranium
        OrbitalParams(92, 0, 0, 7,  -5.0, 2.195659,      0.0,      1.0,      0.0)
        OrbitalParams(92, 0, 1, 7,  -3.0, 1.754808,      0.0,      1.0,      0.0)
        OrbitalParams(92, 0, 2, 6, -10.0, 0.474941, 2.369037, 2.017543, 0.628875)
        OrbitalParams(92, 0, 3, 5,  -7.0, 0.361033, 2.599481, 3.165843, 1.588993);;
        # Neptunium
        OrbitalParams(93, 0, 0, 7,  -5.0, 2.307165,      0.0,      1.0,      0.0)
        OrbitalParams(93, 0, 1, 7,  -3.0, 1.767612,      0.0,      1.0,      0.0)
        OrbitalParams(93, 0, 2, 6, -10.0, 0.928348, 2.504201, 0.106626, 0.132637)
        OrbitalParams(93, 0, 3, 5,  -7.0, 0.278041, 2.636355, 4.100459, 1.036528);;
        # Plutonium
        OrbitalParams(94, 0, 0, 7,  -5.0, 2.418671,      0.0,      1.0,      0.0)
        OrbitalParams(94, 0, 1, 7,  -3.0, 1.780416,      0.0,      1.0,      0.0)
        OrbitalParams(94, 0, 2, 6, -10.0, 1.060298, 2.815564, 0.001859, 0.001426)
        OrbitalParams(94, 0, 3, 5,  -7.0, 0.426227, 2.967294, 2.850708, 0.984260);;
        # Americium
        OrbitalParams(95, 0, 0, 7, -0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(95, 0, 1, 7, -0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(95, 0, 2, 6, -0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(95, 0, 3, 5, -0.0, 0.0, 0.0, 1.0, 0.0);;
        # Curium
        OrbitalParams(96, 0, 0, 7, -0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(96, 0, 1, 7, -0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(96, 0, 2, 6, -0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(96, 0, 3, 5, -0.0, 0.0, 0.0, 1.0, 0.0);;
        # Berkelium
        OrbitalParams(97, 0, 0, 7, -0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(97, 0, 1, 7, -0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(97, 0, 2, 6, -0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(97, 0, 3, 5, -0.0, 0.0, 0.0, 1.0, 0.0);;
        # Californium
        OrbitalParams(98, 0, 0, 7, -0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(98, 0, 1, 7, -0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(98, 0, 2, 6, -0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(98, 0, 3, 5, -0.0, 0.0, 0.0, 1.0, 0.0);;
        # Einsteinium
        OrbitalParams(99, 0, 0, 7, -0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(99, 0, 1, 7, -0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(99, 0, 2, 6, -0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(99, 0, 3, 5, -0.0, 0.0, 0.0, 1.0, 0.0);;
        # Fermium
        OrbitalParams(100, 0, 0, 7, -0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(100, 0, 1, 7, -0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(100, 0, 2, 6, -0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(100, 0, 3, 5, -0.0, 0.0, 0.0, 1.0, 0.0);;
        # Mendeleevium
        OrbitalParams(101, 0, 0, 7, -0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(101, 0, 1, 7, -0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(101, 0, 2, 6, -0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(101, 0, 3, 5, -0.0, 0.0, 0.0, 1.0, 0.0);;
        # Nobelium
        OrbitalParams(102, 0, 0, 7, -0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(102, 0, 1, 7, -0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(102, 0, 2, 6, -0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(102, 0, 3, 5, -0.0, 0.0, 0.0, 1.0, 0.0);;
        # Lawrencium
        OrbitalParams(103, 0, 0, 7, -0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(103, 0, 1, 7, -0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(103, 0, 2, 6, -0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(103, 0, 3, 5, -0.0, 0.0, 0.0, 1.0, 0.0);;
        # Rutherfordium
        OrbitalParams(104, 0, 0, 7, -0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(104, 0, 1, 7, -0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(104, 0, 2, 6, -0.0, 0.0, 0.0, 1.0, 0.0)
        OrbitalParams(104, 0, 3, 5, -0.0, 0.0, 0.0, 1.0, 0.0);;
    ]
)
