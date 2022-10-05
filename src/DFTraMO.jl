module DFTraMO

using LinearAlgebra
using StaticArrays
using Xtal

include("Psphere.jl")
export writePsphere

include("inputs.jl")
export create_run, extract_VASP

end