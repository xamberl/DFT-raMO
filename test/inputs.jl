@testset "Inputs" begin
    cd("input/") do
        ramoinput = read_yaml("ScAl3.yaml")
        @test length(ramoinput.runlist) == 5
        @test ramoinput.emin == -100*DFTraMO.Electrum.EV2HARTREE
        @test ramoinput.emax == 0.25
    end

    @test DFTraMO.parse_sites(["10", "10", "10:2:16"]) == [10, 12, 14, 16]
end
