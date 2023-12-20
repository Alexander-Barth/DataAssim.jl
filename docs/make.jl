using Documenter
using DataAssim

CI = get(ENV, "CI", nothing) == "true"

makedocs(
    format = Documenter.HTML(),
    modules = [DataAssim],
    sitename = "DataAssim",
    pages = [
        "index.md"],
    checkdocs=:none,
)


if CI
    deploydocs(
        repo = "github.com/Alexander-Barth/DataAssim.jl.git",
        target = "build",
    )
end
