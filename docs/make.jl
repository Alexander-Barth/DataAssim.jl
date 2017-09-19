using Documenter
using DataAssim

makedocs(
    format = :html,
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
    julia  = "0.6",
    deps = nothing,
    make = nothing,
)
