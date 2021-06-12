module Path_Test

using Xtals
using Test

@testset "Path Tests" begin
    @test rc[:atomic_masses][:He] == 4.0026
    oldpath = rc[:paths][:crystals]
    rc[:paths][:crystals] = joinpath(pwd(), "other_data", "other_crystals")
    Crystal("other_SBMOF-1.cif")
    rc[:paths][:crystals] = oldpath
    @test true # if made it this far :)
end
end
