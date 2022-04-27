
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
@assert VERSION.minor ≥ 6

using Test, Documenter, Xtals, Graphs, MetaGraphs
Xtals.banner()

for testfile ∈ testfiles
    @info "Running test/$testfile"
    @time include(testfile)
end

# run doctests unless disabled via environment variable
if "doctest" ∉ keys(ENV) || ENV["doctest"] ≠ "false"
    @time doctest(Xtals)
end

@info "Done."
