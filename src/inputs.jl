"""
    create_run(run_name::AbstractString)

Generates a generic input file ("run_name.in") loaded with options for the run in the current directory.
"""
function create_run(run_name::AbstractString)
    open(string(run_name,".in"),"w") do io
    println(io,
        "run_type 1\nsites 1:36\nsite_list sitelist.txt\nAO 1\nradius 2.2\ncontinue_from lastrun.mat")
    end
end

"""
extract_VASP

Searches in the current directory for VASP files and extracts relevant information.
"""
function extract_VASP()
    # Read Fermi energy from OUTCAR
    fermi = get_fermi("OUTCAR")
    
    # Get geometry from POSCAR and KPOINTS (WIP)
    geo = readPOSCAR("POSCAR")
    kpt = readKPOINTS("KPOINTS")
    # Creates an AtomList{3} supercell from POSCAR and KPOINTS
    super = supercell(geo.gen,kpt.grid)

    # Load in default number of atomic orbitals, number of electrons [WIP]

    # Prints info into a system.in file
    #==open("system.in","w") do io
        println(io,"fermi_energy ", fermi, "\nenergy_range -100:"," fermi")
    end==#
    return fermi, geo, kpt, super
end