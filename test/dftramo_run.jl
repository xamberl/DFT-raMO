@testset "dftramo_run" begin
    cd("input/") do
        @test dftramo_run("ScAl3_short.yaml") == nothing
        rm("1_Sc-Sc", force=true, recursive=true)

        @test dftramo_run("displaced_aos.yml") == nothing
        rm("Sc_custom_dx2y2", force=true, recursive=true)
    end
end
