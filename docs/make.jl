#! /usr/bin/env julia
using Documenter
using DFTraMO

makedocs(
    sitename = "DFTraMO",
    format = Documenter.HTML(),
    modules = [DFTraMO]
)

makedocs(
    sitename = "DFT-raMO",
    format = Documenter.HTML(prettyurls = (get(ENV, "CI", nothing) == true)),
    modules = [DFTraMO],
    pages = [
        "Home" => "index.md",
        "API" => "api/index.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#