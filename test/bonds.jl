module Bonds_test

using Xtals, LightGraphs, Test, MetaGraphs, LinearAlgebra

if ! isdir("temp")
    mkdir("temp")
end


function visual_check(xtal::String)
    c = Crystal(xtal)
    c = replicate(c, (2, 2, 2))
    strip_numbers_from_atom_labels!(c)
    infer_geometry_based_bonds!(c, true) # must
    write_xyz(c, "temp/c.xyz")
    write_bond_information(c, "temp/$(c.name)")
    @info c.name * " see .vtk and .xyz to visually check bonds"
end


@testset "bond vectors" begin
    box = replicate(unit_cube(), (10, 10, 10))

    xtal = Crystal(
        "vector test 1",
        box,
        Atoms(
            [:H, :H, :H],
            Frac([
                0.000 0.000 0.075;
                0.075 0.000 0.000;
                0.000 0.000 0.000
            ])
        ),
        Charges{Frac}(0)
    )
    # check non-fatal error for calculating vectors w/o bonds
    @test isnothing(Xtals.calculate_bond_vectors!(xtal))
    infer_bonds!(xtal, true, calculate_vectors=true)
    # check fatal error for attempting to re-calculate vectors
    @test_throws ErrorException calculate_bond_vectors!(xtal)
    # test clearing and re-calculating vectors
    clear_vectors!(xtal)
    @test !get_prop(xtal.bonds, :has_vectors)
    calculate_bond_vectors!(xtal)
    @test get_prop(xtal.bonds, :has_vectors)
    # check vectors
    @test isapprox(get_bond_vector(xtal, 1, 2), [0.00; -0.75;  0.00])
    @test isapprox(get_bond_vector(xtal, 2, 3), [0.75;  0.00;  0.00])
    # check reversal of direction
    @test isapprox(get_bond_vector(xtal, 2, 1), [0.00;  0.75;  0.00])
    # check bond angle
    @test isapprox(bond_angle(xtal, 1, 2, 3), deg2rad(90))
    # check invalid bond angle
    @test_throws AssertionError bond_angle(xtal, 2, 3, 1)
    # check bond distance calculations
    @test isapprox(bond_distance(xtal, 1, 2), norm(get_prop(xtal.bonds, 1, 2, :vector)))

    xtal = Crystal(
        "vector test 2",
        box,
        Atoms(
            [:H, :H, :H],
            Frac([
                1.000 0.000 0.075;
                0.075 0.000 0.000;
                0.000 0.000 0.000
            ])
        ),
        Charges{Frac}(0)
    )
    infer_bonds!(xtal, true, calculate_vectors=true)
    @test isapprox(get_bond_vector(xtal, 1, 2), [0.00; -0.75;  0.00])
    @test isapprox(get_bond_vector(xtal, 2, 1), [0.00;  0.75;  0.00])
    @test isapprox(get_bond_vector(xtal, 2, 3), [0.75;  0.00;  0.00])
    @test isapprox(bond_angle(xtal, 1, 2, 3), deg2rad(90)) # 90° angles work to numerical precision

    xtal = Crystal(
        "vector test 3",
        box,
        Atoms(
            [:H, :H, :H],
            Frac([
                0.000 0.075 0.150;
                0.000 0.000 0.000;
                0.000 0.000 0.000
            ])
        ),
        Charges{Frac}(0)
    )
    infer_bonds!(xtal, true, calculate_vectors=true)
    @test isapprox(get_bond_vector(xtal, 1, 2), [0.75; 0.00; 0.00])
    @test isapprox(get_bond_vector(xtal, 2, 3), [0.75; 0.00; 0.00])
    @test isapprox(bond_angle(xtal, 1, 2, 3), π, rtol=1e-4) # 180° angles work to within 0.01%

    methane = Frac(read_xyz(joinpath(rc[:paths][:crystals], "methane.xyz")), box)
    xtal = Crystal(
        "vector test 4",
        box,
        methane,
        Charges{Frac}(0)
    )
    infer_bonds!(xtal, true, calculate_vectors=true)
    # test that ∠ijk = ∠kji and all angles are all the same
    α = bond_angle(xtal, 2, 1, 3)
    results = Bool[]
    for i ∈ 2:5
        for j ∈ 2:5
            if i == j
                push!(results, isapprox(bond_angle(xtal, i, 1, j), 0)) # 0° angles work to numerical precision
            else
                push!(results, isapprox(bond_angle(xtal, i, 1, j), α, rtol=1e-5)) # identical angles work to within 0.001%
            end
        end        
    end
    @test all(results)
    # test that bond angles match expected for tetrahedral carbon
    @test isapprox(bond_angle(xtal, 2, 1, 3), deg2rad(109.5), rtol=1e-3) # 109.5° angles work to within 0.1%
    # test that result is translation invariant
    results = Vector{Bool}(undef, 10000)
    for i ∈ 1:length(results)
        translate_by!(xtal.atoms.coords, Frac(rand(3, 1)))
        results[i] = isapprox(bond_angle(xtal, 2, 1, 3), α, rtol=1e-5)
    end
    @test all(results)
end
@testset "bond sanity check" begin
    box = unit_cube()
    texas_carbon = read_xyz(joinpath(rc[:paths][:crystals], "texas_carbon.xyz"))
    xtal = Crystal("CH5", box, Atoms(texas_carbon.n, texas_carbon.species, Frac(texas_carbon.coords, box)), Charges{Frac}(0))
    @test !bond_sanity_check(xtal)[1]
    infer_bonds!(xtal, false)
    @test !bond_sanity_check(xtal)[1]
    trihydrogen = read_xyz(joinpath(rc[:paths][:crystals], "trihydrogen.xyz"))
    xtal = Crystal("H3", box, Atoms(trihydrogen.n, trihydrogen.species, Frac(trihydrogen.coords, box)), Charges{Frac}(0))
    infer_bonds!(xtal, false)
    @test !bond_sanity_check(xtal)[1]
end
@testset "NiPyC2 Tests" begin
    # NiPyC2 bonding
    c = Crystal("NiPyC2_relax.cif")
    strip_numbers_from_atom_labels!(c)
    c = replicate(c, (4, 4, 4))
    bonding_rules = [BondingRule(:H, :*, 1.2),
                     BondingRule(:Ni, :O, 2.5),
                     BondingRule(:Ni, :N, 2.5),
                     BondingRule(:*, :*, 1.9)]
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
    covalent_radii = rc[:covalent_radii]
    covalent_radii[:Cu] = 1.15
    remove_bonds!(c)

    @test ne(c.bonds) == 0

    infer_geometry_based_bonds!(c, true, covalent_radii=covalent_radii)

    @test c.atoms.species[neighbors(c.bonds, 1)] == [:O, :O, :O, :O]
end
@testset "BACMOH Tests" begin
    xtal = Crystal("BACMOH_clean.cif")
    strip_numbers_from_atom_labels!(xtal)
    infer_bonds!(xtal, true)

    @test bond_sanity_check(xtal)
    @test xtal.atoms.species[1] == :Cu
    @test xtal.atoms.species[neighbors(xtal.bonds, 1)] == [:O, :N, :N, :O]
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
    scipy = rc[:scipy]
    rc[:scipy] = nothing
    @test_throws ErrorException infer_geometry_based_bonds!(Crystal("IRMOF-1.cif"), true)
    rc[:scipy] = scipy
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
@testset "bonding rules" begin
    write_bonding_rules("temp/bonding_rules.csv")
    @test isfile("temp/bonding_rules.csv")
    bonding_rules = read_bonding_rules("temp/bonding_rules.csv")
    n = 4656 # the number generated by default. subject to change.
    @test length(bonding_rules) == n
    add_bonding_rules([BondingRule(:foo, :bar, 5.)])
    @test length(rc[:bonding_rules]) == (n + 1)
    println([r for r ∈ bonding_rules if :C == r.species_j][1])
    @test true
    rc[:bonding_rules] = [BondingRule(:foo, :bar, 0.0)]
    xtal = Crystal("cof-102.cif")
    @test isnan(Xtals.is_bonded(xtal, 1, 2, rc[:bonding_rules])[2])
    rc[:bonding_rules] = bondingrules()
    @test true
end
@testset "etc" begin
    xtal = Crystal("SBMOF-1.cif")
    write_bond_information(xtal, "temp/nothing.vtk", center_at_origin=true)
    @test isfile("temp/nothing.vtk")
    xtal = Crystal("IRMOF-1.cif")
    infer_bonds!(xtal, true)
    write_xyz(xtal, "temp/IRMOF-1.xyz")
    write_bond_information(xtal, "temp/all_bonds.vtk")
    write_bond_information(xtal, "temp/no_pb.vtk", bond_filter=:cross_boundary=>p->!p)
end
@testset "bonds from xyz" begin
    xtal = Crystal("IRMOF-1.cif")
    infer_bonds!(xtal, false)
    bonds = infer_bonds(Cart(xtal.atoms, xtal.box))
    @test bonds == xtal.bonds
end
end
# visual_check
