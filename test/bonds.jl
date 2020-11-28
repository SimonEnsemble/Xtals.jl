module Bonds_test
using Xtals
using LightGraphs
using Test

function visual_check(xtal::String)
    c = Crystal(xtal)
    c = replicate(c, (2, 2, 2))
    strip_numbers_from_atom_labels!(c)
    infer_geometry_based_bonds!(c, true) # must
    write_xyz(c)
    write_bond_information(c)
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
    infer_bonds!(c, true, bonding_rules)
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
    write_bond_information.([c, c1], ["c", "c1"])
    write_xyz(c, "c")

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
    cordero_params = cordero_parameters()
    cordero_params[:Cu] = Dict(:radius_Å => 1.15, :esd_pm => 4.)
    remove_bonds!(c)

    @test ne(c.bonds) == 0

    infer_geometry_based_bonds!(c, true, cordero_params=cordero_params)

    @test c.atoms.species[neighbors(c.bonds, 1)] == [:O, :O, :O, :O]
end
end
