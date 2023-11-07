@testset OccupiedStates begin
    wf = readWAVECAR("input/WAVECAR")
    # Fermi energy is the maximum...obtain it separately
    emax = get_fermi("input/OUTCAR")[:fermi]
end
