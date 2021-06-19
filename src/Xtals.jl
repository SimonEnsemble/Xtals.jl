module Xtals

using Bio3DView, CSV, DataFrames, LightGraphs, LinearAlgebra, MetaGraphs, Printf, PyCall, UUIDs

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
include("cpk_colors.jl")
include("covalent_radii.jl")
include("bonds.jl")

rc[:atomic_masses] = DEFAULT_ATOMIC_MASSES
rc[:cpk_colors] = DEFAULT_CPK_COLORS
rc[:covalent_radii] = DEFAULT_COVALENT_RADII
DEFAULT_BONDING_RULES = bondingrules()
rc[:bonding_rules] = DEFAULT_BONDING_RULES

# - lists python dependencies
# - checks for python dependencies at pre-compile, warns on failure
# - provides function for python dependency instantiation
include("pydeps.jl")


function __init__()
    # load Python dependencies
    init_pydeps()
    # create path entries in global dictionary
    rc[:paths][:crystals] = ""
    rc[:paths][:data] = ""
    # set paths to data and crystals relative to pwd() at import
    set_paths(joinpath(pwd(), "data"))
end


export
    # Xtals.jl
    rc,

    # matter.jl
    Coords, Frac, Cart, Atoms, Charges, wrap!, neutral, net_charge, translate_by!,

    # box.jl
    Box, replicate, unit_cube, write_vtk, inside, fractional_coords, cartesian_coords,

    # distance.jl
    nearest_image!, distance, overlap, remove_duplicates, pairwise_distances,

    # misc.jl
    read_xyz, write_xyz, read_mol, write_mol2, assert_P1_symmetry, set_paths, view_crystal,

    # crystal.jl
    Crystal, strip_numbers_from_atom_labels!, assign_charges, chemical_formula, molecular_weight, 
    crystal_density, write_cif, has_charges, apply_symmetry_operations, write_cssr, primitive_cell,

    # bonds.jl
    infer_bonds!, write_bond_information, BondingRule, bond_sanity_check, remove_bonds!, 
    infer_geometry_based_bonds!, get_bonding_rules, set_bonding_rules, read_bonding_rules,
    write_bonding_rules, add_bonding_rules, drop_cross_pb_bonds!, bondingrules, bond_angle,
    get_bond_vector, calculate_bond_vectors!, clear_vectors!, bond_distance

end # module Xtals
