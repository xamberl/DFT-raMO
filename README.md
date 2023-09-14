# DFT-raMO

[![Documentation (stable)][docs-stable-img]][docs-stable-url]
[![Documentation (dev)][docs-dev-img]][docs-dev-url]
[![CI status][ci-status-img]][ci-status-url]

[![DFTraMO.jl logo][dftramo-logo]][dftramo-url]

DFT-raMO, a Wannier-type analysis, and its supporting tools.

# The original DFT-raMO (written in MATLAB)

The current version of DFT-raMO is a MATLAB script which can be run in the terminal. You'll need
a MATLAB installation in order to run it.

# DFTraMO.jl
A Julia package that performs the same DFT-raMO analysis, but with support for parallelization.
The results from this package are higher in quality than the original MATLAB code, and it has
support for both abinit and VASP inputs. Runs are specified through YAML files
([an example](examples/ScAl3.yaml) is provided).

This package also contains Psphere analysis, a tool to determine the degree of localization for
reconstructed orbitals.

## Installation
This package is not currently registered in the Julia package registry, so you'll need to add the
package by the URL directly:
```
(v1.8+ pkg)> add https://github.com/xamberl/DFT-raMO#main
```
The `#main` portion is optional, and can be changed to track a different branch.

[docs-stable-img]:  https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]:  https://xamberl.github.io/DFT-raMO/stable
[docs-dev-img]:     https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]:     https://xamberl.github.io/DFT-raMO/dev
[ci-status-img]:    https://github.com/xamberl/DFT-raMO/workflows/CI/badge.svg
[ci-status-url]:    https://github.com/xamberl/DFT-raMO/actions
[dftramo-logo]:     docs/src/assets/dftramologo_500.png
[dftramo-url]:      https://github.com/xamberl/DFT-raMO
