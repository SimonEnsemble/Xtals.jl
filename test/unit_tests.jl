testfiles = [
    "bonds.jl",
    "misc.jl",
    "crystal.jl",
    "matter.jl",
    "distance.jl",
    "box.jl",
    "assert_p1_symmetry.jl",
    "paths.jl"
    ]

using Test, Graphs, MetaGraphs, Documenter

if !isdir("temp")
    mkdir("temp")
end

using Xtals

Xtals.banner()

for testfile ∈ testfiles
    @info "Running test/$testfile"
    include(testfile)
end

@info "Tests complete!"
