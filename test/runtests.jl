testfiles = [
    "matter.jl",
    "distance.jl",
    "crystal.jl",
    "misc.jl",
    "bonds.jl",
    "box.jl",
    "assert_p1_symmetry.jl",
    "paths.jl"
    ]

using Test, LightGraphs, MetaGraphs
using Xtals

@info "\n\n\t\tXtals.jl\n\n\n"

if !isdir("temp")
    mkdir("temp")
end

include.(testfiles)

@info "Done!"
