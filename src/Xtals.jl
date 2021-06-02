module Xtals

using CSV, DataFrames, Printf, LinearAlgebra, LightGraphs, PyCall, MetaGraphs

# global variable dictionary
global rc = Dict{Symbol,Any}()
rc[:paths] = Dict{Symbol,String}()

include("matter.jl")
include("box.jl")
include("crystal.jl")
include("distance.jl")
include("misc.jl")
include("repfactors.jl")

include("atomic_masses.jl")
rc[:atomic_masses] = DEFAULT_ATOMIC_MASSES
include("cpk_colors.jl")
rc[:cpk_colors] = DEFAULT_CPK_COLORS
include("covalent_radii.jl")
rc[:covalent_radii] = DEFAULT_COVALENT_RADII
include("bonds.jl")
rc[:bonding_rules] = DEFAULT_BONDING_RULES


function load_pydep(pydep)
    try
        return pyimport(pydep)
    catch
        @warn "Error loading $pydep. Some functionaltiy may be missing."
        return nothing
    end
end


function __init__()
    # load Python dependencies
    rc[:scipy] = load_pydep("scipy.spatial")
    
    # sets paths to data and crystals relative to pwd() at import
    rc[:paths][:data] = joinpath(pwd(), "data")
    rc[:paths][:crystals] = joinpath(rc[:paths][:data], "crystals")
    @debug "Environment variables" rc
end


export
    # Xtals.jl
    rc,

    # matter.jl
    Coords, Frac, Cart, Atoms, Charges, wrap!, neutral, net_charge,
    translate_by!, origin,

    # box.jl
    Box, replicate, unit_cube, write_vtk, inside, fractional_coords,
    cartesian_coords,

    # distance.jl
    nearest_image!, distance, overlap, remove_duplicates, pairwise_distances,

    # misc.jl
    read_xyz, write_xyz, read_mol, write_mol2, assert_P1_symmetry,

    # crystal.jl
    Crystal, strip_numbers_from_atom_labels!, assign_charges,
    chemical_formula, molecular_weight, crystal_density, write_cif, has_charges,
    apply_symmetry_operations, write_cssr,

    # bonds.jl
    infer_bonds!, write_bond_information, BondingRule, bond_sanity_check,
    remove_bonds!, infer_geometry_based_bonds!, get_bonding_rules,
    set_bonding_rules, read_bonding_rules, write_bonding_rules, add_bonding_rules,
    drop_cross_pb_bonds!, bondingrules

end # module Xtals
