@testset "OccupiedStates" begin
    wf = readWAVECAR("input/WAVECAR")
    # Fermi energy is the maximum...obtain it separately
    emax = get_fermi("input/OUTCAR")[:fermi]
    occ_states = OccupiedStates(wf, emin = min_energy(wf), emax = emax)
    # There shouldn't be more G-vectors than the wavefunction bounds permit
    @test size(occ_states, 1) <= prod(size(wf)[4:end])
    @test size(occ_states, 1) == length(occ_states.G)
    for g in occ_states.g
        @test g in SVector.(Tuple.(FFTBins(wf)))
    end
    # Size of the second dimension depends on the energy range
    @test size(occ_states, 2) == count(emin .<= wf.energies .<= emax)
    @test size(occ_states, 2) == length(occ_states.skb)
    # SpinKPointBand structs should return all relevant info from wavefunction
    for x in occ_states.skb
        @test x.spin in wf.spins
        @test KPoint(x) in KPointMesh(wf)
        @test x.band in axes(wf, 3)
    end
    # Not sure if encoding equality is worth doing
    @test occ_states == OccupiedStates(wf, emax = emax)
end
