testfiles = [
    "misc.jl",
    "bonds.jl",
    "matter.jl",
    "distance.jl",
    "crystal.jl",
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

for testfile ∈ testfiles
    @info "Running test/$testfile"
    include(testfile)
end

@info "Done!"
