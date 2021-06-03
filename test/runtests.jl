testfiles = [
    "bonds.jl",
    "matter.jl",
    "distance.jl",
    "crystal.jl",
    "misc.jl",
    "box.jl",
    "assert_p1_symmetry.jl",
    "paths.jl"
    ]

using Test, LightGraphs, MetaGraphs, Documenter

using Xtals

@info "\n\n\t\tXtals.jl\n\n\n"

@info "Running unit tests..."
for testfile ∈ testfiles
    @info "Running test/$testfile"
    include(testfile)
end

# run doctests
doctest(Xtals)

@info "Done!"
