@testset "dftramo_run" begin
    cd("input/") do
        @test dftramo_run("ScAl3_short.yaml") == nothing
        rm("1_Sc-Sc", force=true, recursive=true)
    end
end
