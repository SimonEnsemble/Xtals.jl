version_pattern = "[0-9].[0-9].[0-9]"
binary_path_string = match(Regex("julia-$version_pattern/bin/julia"), ENV["_"]).match
version_string = match(Regex("$version_pattern"), binary_path_string).match
major_version, minor_version, _ = parse.(Int, split(version_string, "."))

@assert major_version == 1
@assert minor_version ≥ 5
include("unit_tests.jl")


if minor_version == 7
    include("doc_tests.jl")
end
