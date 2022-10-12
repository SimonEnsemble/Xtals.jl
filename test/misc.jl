module Misc_Test

using Xtals, Test, DataFrames, Graphs, CSV
import IOCapture.capture

@testset "Misc Tests" begin
    xtal = Crystal("SBMOF-1.cif")
    infer_bonds!(xtal, true)
    capture() do
        @test_throws String view_crystal(xtal) # this function really can't be tested programmatically...
    end

    am = rc[:atomic_masses]
    @test isapprox(am[:H], 1.00794, atol=0.001)
    @test isapprox(am[:Co], 58.9332, atol=0.001)

    test_xyz_filename = "atoms_test"
    c = Cart([
        1.0 4.0
        2.0 5.0
        3.0 6.0
    ])
    s = [:C, :H]
    atoms = Atoms(s, c)
    write_xyz(atoms, test_xyz_filename)
    atoms_read = read_xyz(test_xyz_filename * ".xyz")
    @test isapprox(atoms, atoms_read)
    rm(test_xyz_filename * ".xyz")

    @test rc[:cpk_colors][:Li] == (204, 128, 255)

    atoms, bonds, bond_types = read_mol("data/example.mol")
    @test (atoms.species[1] == :O) &&
          (atoms.species[end] == :H) & (length(atoms.species) == 31) &&
          (atoms.n == 31)
    @test ne(bonds) == 32
    @test nv(bonds) == 31
    @test (bond_types[1] == 1) && (bond_types[3] == 2)
    @test neighbors(bonds, 1) == [2, 5]

    xtal = Crystal("SBMOF-1.cif"; infer_bonds=true, periodic_boundaries=true)
    mol2_temp = tempname()
    write_mol2(xtal; filename=mol2_temp)
    @test isfile(mol2_temp)
    write_mol2(xtal)
    @test isfile("SBMOF-1.mol2")
    write_xyz(xtal)
    @test isfile("SBMOF-1.xyz")
    write_xyz(xtal, "foo.xyz")
    @test isfile("foo.xyz")
end

end
