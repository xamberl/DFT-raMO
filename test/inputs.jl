@testset "Inputs" begin
    cd("input/") do
        ramoinput = read_yaml("ScAl3.yaml")
        @test length(ramoinput.runlist) == 5
        @test ramoinput.emin == -100*DFTraMO.Electrum.EV2HARTREE
        @test ramoinput.emax == 0.25
        @test ramoinput.mode == "auto_psphere"
        @test ramoinput.runlist[5].radius == 3.97*DFTraMO.Electrum.BOHR2ANG

        ramoinput = read_yaml("displaced_aos.yml")
        @test ramoinput.mode == "discard"
        @test ramoinput[1].sites == [1, 3, 7, 8]
        @test ramoinput[1].type == "dx2y2"
        @test ramoinput[1].direction == [0.5, 0.5, 0]/DFTraMO.norm([0.5, 0.5, 0])
        @test ramoinput[1].radius == 0.5
    end

    @test DFTraMO.parse_sites(["10", "10", "10:2:16"]) == [10, 12, 14, 16]
end
