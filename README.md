# DFT-raMO
DFT-raMO and tools to support DFT-raMO, a Wannier-type analysis.

# The original DFT-raMO (written in MATLAB)

The current version of DFT-raMO is a MATLAB script which can be run in the terminal. You'll need
a MATLAB installation in order to run it.

## Psphere.jl
A Julia script to batch-process raMOs from DFT-raMO in order to obtain Psphere, a metric describing
the degree of localization of a raMO function.

### Usage
You will need the LinearAlgebra package. You can add this in the Julia package manager.
In the Julia REPL, run the Psphere.jl code to export the function with `include("Psphere.jl")`.
Now you may use the `writePsphere` function in the REPL. Type `?writePsphere` to see inputs.
