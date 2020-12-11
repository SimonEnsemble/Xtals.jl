module Misc_Test

using Xtals
using Test
using DataFrames
using LightGraphs
using CSV

@testset "Misc Tests" begin
    am = get_atomic_masses()
    @test isapprox(am[:H], 1.00794, atol=0.001)
    @test isapprox(am[:Co], 58.9332, atol=0.001)

    test_xyz_filename = "atoms_test"
    c = Cart([1.0 4.0;
              2.0 5.0;
              3.0 6.0]
             )
    s = [:C, :H]
    atoms = Atoms(s, c)
    write_xyz(atoms, test_xyz_filename)
    atoms_read = read_xyz(test_xyz_filename * ".xyz")
    @test isapprox(atoms, atoms_read)
    rm(test_xyz_filename * ".xyz")

    @test get_cpk_colors()[:Li] == (204,128,255)

    atoms, bonds, bond_types = read_mol("example.mol")
    @test (atoms.species[1] == :O) && (atoms.species[end] == :H) & (length(atoms.species) == 31) && (atoms.n == 31)
    @test ne(bonds) == 32
    @test nv(bonds) == 31
    @test (bond_types[1] == 1) && (bond_types[3] == 2)
    @test neighbors(bonds, 1) == [2, 5]
end
end
