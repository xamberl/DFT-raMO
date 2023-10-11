@testset "Inputs" begin
    cd("input/") do
        (runs, checkpoint, auto_psphere, dftinfo, emin, emax) = read_run_yaml("ScAl3.yaml")
        @test length(runs) == 5
        @test emin == -100*DFTraMO.Electrum.EV2HARTREE
        @test emax == 0.25
    end

    @test DFTraMO.parse_sites(["1:10"]) == collect(1:10)
    @test DFTraMO.parse_sites(["10", "10", "10:2:16"]) == [10, 12, 14, 16]
end
