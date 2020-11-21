module Path_Test

using Xtals
using Test

@testset "Path Tests" begin
    # recommended way to change path to files
    @eval Xtals PATH_TO_CRYSTALS = joinpath(pwd(), "other_data", "other_crystals")
    Crystal("other_SBMOF-1.cif")
    @test true # if made it this far :)
end
end
