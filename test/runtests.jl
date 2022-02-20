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

@assert VERSION.major == 1
@assert VERSION.minor ≥ 1

using Test, Documenter, Xtals, Graphs, MetaGraphs
Xtals.banner()

if !isdir("temp")
    mkdir("temp")
end

for testfile ∈ testfiles
    @info "Running test/$testfile"
    @time include(testfile)
end

if VERSION.minor ≥ 7
    @time doctest(Xtals)
end

@info "Done."
