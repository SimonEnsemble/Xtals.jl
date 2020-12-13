module Bonds_test
using Xtals, LightGraphs, Test, MetaGraphs

function visual_check(xtal::String)
    c = Crystal(xtal)
    c = replicate(c, (2, 2, 2))
    strip_numbers_from_atom_labels!(c)
    infer_geometry_based_bonds!(c, true) # must
    write_xyz(c, "temp/c.xyz")
    write_bond_information(c, "temp/$(c.name)")
    @info c.name * " see .vtk and .xyz to visually check bonds"
end

@testset "NiPyC2 Tests" begin
    # NiPyC2 bonding
    c = Crystal("NiPyC2_relax.cif")
    strip_numbers_from_atom_labels!(c)
    c = replicate(c, (4, 4, 4))
    bonding_rules = [BondingRule(:H, :*, 0.4, 1.2),
                     BondingRule(:Ni, :O, 0.4, 2.5),
                     BondingRule(:Ni, :N, 0.4, 2.5),
                     BondingRule(:*, :*, 0.4, 1.9)]
    infer_bonds!(c, true, bonding_rules=bonding_rules)
    c1 = deepcopy(c)
    conn_comps = connected_components(c.bonds)

    @test length(conn_comps) == 2 # interpenetrated

    c_red = getindex(c, conn_comps[1])
    c_blue = getindex(c, conn_comps[2])

    @test ne(c_red.bonds) + ne(c_blue.bonds) == ne(c.bonds)

    remove_bonds!(c)
    infer_geometry_based_bonds!(c, true)

    @test length(conn_comps) == 2 # interpenetrated

    # debugging outputs (delete me)
    write_bond_information.([c, c1], ["temp/c", "temp/c1"])
    write_xyz(c, "temp/c")

    @test c1.bonds == c.bonds # consistency between two different bonding schemes

end
@testset "FIQCEN Tests" begin
    # FIQCEN bonding
    c = Crystal("FIQCEN_clean.cif")
    strip_numbers_from_atom_labels!(c)
    infer_geometry_based_bonds!(c, true)

    @test length(connected_components(c.bonds)) == 1 # not interpenetrated

    @test c.atoms.species[neighbors(c.bonds, 1)] == [:Cu, :O, :O, :O, :O]

    visual_check("FIQCEN_clean.cif")

    # reduce covalant radius to see Cu-Cu bond disappear
    covalent_radii = Xtals.get_covalent_radii()
    covalent_radii[:Cu] = Dict(:radius_Å => 1.15, :esd_pm => 4.)
    remove_bonds!(c)

    @test ne(c.bonds) == 0

    infer_geometry_based_bonds!(c, true, covalent_radii=covalent_radii)

    @test c.atoms.species[neighbors(c.bonds, 1)] == [:O, :O, :O, :O]
end
@testset "bond inference options" begin
    xtal1 = Crystal("SBMOF-1_overlap.cif", remove_duplicates=true)
    # infer_bonds option tests
    xtal2 = deepcopy(xtal1) # no bonds
    infer_bonds!(xtal1, true)
    infer_bonds!(xtal2, false)
    xtal3 = Crystal("SBMOF-1_overlap.cif", remove_duplicates=true,
    infer_bonds=:cordero, periodic_boundaries=true)
    xtal4 = Crystal("SBMOF-1_overlap.cif", remove_duplicates=true,
    infer_bonds=:cordero, periodic_boundaries=false)
    xtal5 = Crystal("SBMOF-1_overlap.cif", remove_duplicates=true,
    infer_bonds=:voronoi, periodic_boundaries=true)
    xtal6 = Crystal("SBMOF-1_overlap.cif", remove_duplicates=true,
    infer_bonds=:voronoi, periodic_boundaries=false)
    # xtal1  : default bonding rules, periodic boundaries
    # xtal2 : default bonding rules, no periodic boundaries
    # xtal3 : default bonding rules, periodic boundaries
    # xtal4 : default bonding rules, no periodic boundaries
    # xtal5 : voronoi bonding rules, periodic boundaries
    # xtal6 : voronoi bonding rules, no periodic boundaries
    @test xtal1.bonds ≠ xtal2.bonds
    @test xtal1.bonds == xtal3.bonds
    @test xtal1.bonds ≠ xtal4.bonds
    # @test xtal1.bonds == xtal5.bonds # TODO examine test failure
    @test xtal1.bonds ≠ xtal6.bonds
    @test xtal2.bonds ≠ xtal3.bonds
    @test xtal2.bonds == xtal4.bonds
    @test xtal2.bonds ≠ xtal5.bonds
    #@test xtal2.bonds == xtal6.bonds # TODO examine test failure
    @test xtal3.bonds ≠ xtal4.bonds
    #@test xtal3.bonds == xtal5.bonds # TODO examine test failure
    @test xtal3.bonds ≠ xtal6.bonds
    @test xtal4.bonds ≠ xtal5.bonds
    #@test xtal4.bonds == xtal6.bonds # TODO examine test failure
    @test xtal5.bonds ≠ xtal6.bonds
end
@testset ".mol/.cif bonds vs. inferred" begin
    mol_atoms, mol_bonds = read_mol("data/example.mol")
    box = unit_cube()
    xtal = Crystal("example.mol", box, Frac(mol_atoms, box), Charges{Frac}(0))
    infer_bonds!(xtal, false)
    @test mol_bonds == xtal.bonds
    mol_atoms, mol_bonds = read_mol("data/cof-102.mol")
    xtal = Crystal("cof-102.cif")
    infer_bonds!(xtal, false) # source mol file has no periodic bonds
    @test mol_bonds == xtal.bonds
end
@testset "metadata" begin
    xtal = Crystal("cof-102.cif")
    infer_bonds!(xtal, true)
    bonds = deepcopy(xtal.bonds)
    e = collect(edges(xtal.bonds))[1]
    set_prop!(xtal.bonds, src(e), dst(e), :distance, missing)
    Xtals.calc_missing_bond_distances!(xtal)
    @test xtal.bonds == bonds
    xtal2 = Crystal("cof-102.cif")
    infer_bonds!(xtal2, false)
    drop_cross_pb_bonds!(xtal)
    @test xtal.bonds == xtal2.bonds
    write_mol2(xtal, filename="temp/cof-102_no_pb.mol2")
end
@testset "bonding rules api" begin
    write_bonding_rules("temp/bonding_rules.csv")
    @test isfile("temp/bonding_rules.csv")
    bonding_rules = read_bonding_rules("temp/bonding_rules.csv")
    n = 5565 # the number generated by default. subject to change.
    @test length(bonding_rules) == n
    add_bonding_rules([BondingRule(:foo, :bar, 0., 5.)])
    @test length(get_bonding_rules()) == (n + 1)
    println([r for r ∈ bonding_rules if :C == r.species_j][1])
    @test true
    set_bonding_rules(Xtals.bondingrules()) # for re-run safety
end
@testset "etc" begin
    xtal = Crystal("SBMOF-1.cif")
    write_bond_information(xtal, "temp/nothing.vtk", center_at_origin=true)
    @test isfile("temp/nothing.vtk")
end
end
# visual_check
