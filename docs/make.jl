using Documenter
using Xtals

makedocs(
    root = joinpath(dirname(pathof(Xtals)), "..", "docs"),
    modules = [Xtals],
    sitename = "Xtals.jl",
    clean = true,
    pages = [
            "Xtals" => "index.md",
            "matter" => "matter.md",
            "boxes" => "box.md",
            "crystals" => "crystal.md",
            "computing distances" => "distance.md",
            "bonding" => "bonds.md"
            ],
    format = Documenter.HTML(assets = ["assets/flux.css"])
)

deploydocs(
    repo = "github.com/SimonEnsemble/Xtals.jl.git",
 #     push_preview=false,
 #     deps = Deps.pip("mkdocs", "mkdocs-material", "pymdown-extensions") # These are dependencies for the site, not the package
)
