using Test
using StaticArrays
using Electrum
using DFTraMO

@testset "All tests" begin
    # Add your tests here, ideally as file includes.
    include("inputs.jl")
    include("occupiedstates.jl")
    include("dftramo_run.jl")
end
