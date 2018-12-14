using Documenter
using DataAssim

makedocs(
    format = Documenter.HTML(),
    modules = [DataAssim],
    sitename = "DataAssim",
    pages = [
        "index.md"]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.

deploydocs(
    repo = "github.com/Alexander-Barth/DataAssim.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
