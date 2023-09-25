#! /usr/bin/env julia
using Documenter
using DFTraMO


makedocs(
    sitename = "DFT-raMO.jl",
    format = Documenter.HTML(prettyurls = (get(ENV, "CI", nothing) == true)),
    modules = [DFTraMO],
    pages = [
        "Home"      => "index.md",
        "Theory"    => "theory.md",
        "Usage"     => "usage.md",
        "Tutorials"  => ["tutorials/IrIn3.md"],
        "API"       => "api/index.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/xamberl/DFT-raMO.git"
)
