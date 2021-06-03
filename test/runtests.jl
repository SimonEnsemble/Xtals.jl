testfiles = [
    "crystal.jl",
    "misc.jl",
    "bonds.jl",
    "matter.jl",
    "distance.jl",
    "box.jl",
    "assert_p1_symmetry.jl",
    "paths.jl"
    ]

using Test, LightGraphs, MetaGraphs, Documenter

if !isdir("temp")
    mkdir("temp")
end

@info "\n\n\t\tXtals.jl\n\n\n"

using Xtals

for testfile ∈ testfiles
    @info "Running test/$testfile"
    include(testfile)
end

# run doctests
doctest(Xtals)

@info "Tests complete!"
