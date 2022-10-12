testfiles = [
    "crystal.jl"
    "bonds.jl"
    "misc.jl"
    "matter.jl"
    "distance.jl"
    "box.jl"
    "assert_p1_symmetry.jl"
    "paths.jl"
]

@assert (VERSION.major == 1) && (VERSION.minor ≥ 6) "Minimum Julia version not met."

using Documenter, IOCapture, Logging, Test, Xtals

Xtals.banner()

for testfile in testfiles
    @info "Running test/$testfile"
    with_logger(NullLogger()) do
        return include(testfile)
    end
end

# run doctests unless disabled via environment variable
if "doctest" ∉ keys(ENV) || ENV["doctest"] ≠ "false"
    doctest(Xtals)
end

@info "Done."
