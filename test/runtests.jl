@assert VERSION.major == 1
@assert VERSION.minor ≥ 0

include("unit_tests.jl")

if VERSION.minor == 7
    include("doc_tests.jl")
end
