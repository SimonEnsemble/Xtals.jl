module Path_Test

using Xtals
using Test

@testset "Path Tests" begin
    @test rc[:atomic_masses][:He] == 4.0026
    oldpath = rc[:paths][:data]
    newpath = joinpath(pwd(), "other_data")
    set_paths(newpath, print_paths=true)
    @test rc[:paths][:crystals] == joinpath(newpath, "crystals")
    rc[:paths][:crystals] = joinpath(newpath, "other_crystals")
    xtal = Crystal("SBMOF-1.cif")
    @test xtal.atoms.n == 120
    rc[:paths][:foo] = "bar"
    set_paths(oldpath)
    @test rc[:paths][:foo] == joinpath(oldpath, "foo")
end
end
