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

using Test, LightGraphs, MetaGraphs, Documenter, FIGlet

if !isdir("temp")
    mkdir("temp")
end

FIGlet.render("Xtals.jl", FIGlet.availablefonts()[5])

using Xtals

for testfile ∈ testfiles
    @info "Running test/$testfile"
    include(testfile)
end

# run doctests
doctest(Xtals)

@info "Tests complete!"
