# DFT-raMO
DFT-raMO and tools to support DFT-raMO, a Wannier-type analysis.

# The original DFT-raMO (written in MATLAB)

The current version of DFT-raMO is a MATLAB script which can be run in the terminal. You'll need
a MATLAB installation in order to run it.

# DFTraMO.jl
A Julia package that performs the same DFT-raMO analysis, but with support for parallelization.
This code can be adapted to run 

This package also contains Psphere analysis, a tool to determine the degree of localization for
reconstructed orbitals.

## Usage
This package is not currently registered in the Julia package registry, and is located in a private
Git repo as of this writing. You'll need to clone the package from Git over SSH first:

```bash
[user@host]$ git clone git@github.com:xamberl/DFT-raMO
```

Then you can start a Julia instance, enter package mode (right bracket),  and use `activate 
[clone directory]` to use and modify the package.
