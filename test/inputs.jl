@testset "Inputs" begin
    cd("input/") do
        input = read_run_yaml("ScAl3.yaml")
    end
end
