@testset "Inputs" begin
    cd("input/") do
        input = read_run_yaml("ScAl3.yaml")
    end

    @test DFTraMO.parse_sites(["1:10"]) == collect(1:10)
    @test DFTraMO.parse_sites(["10", "10", "10:2:16"]) == [10, 12, 14, 16]
end
