module Xtals

using CSV, DataFrames, Printf, LinearAlgebra, LightGraphs, PyCall, MetaGraphs

include("paths_and_globals.jl")
include("matter.jl")
include("box.jl")
include("crystal.jl")
include("distance.jl")
include("misc.jl")
include("covalent_radii.jl")
include("bonds.jl")
include("repfactors.jl")
include("atomic_masses.jl")
include("cpk_colors.jl")


# runs every time the module is imported
function __init__()
    # global variable dictionary, accessed via get/set_global API
    global GLOBAL = Dict{Symbol,Any}()
    # if the user changes directory, path_to_data will change as well
    set_global(:path_to_data => joinpath(pwd(), "data"))
    set_global(:path_to_crystals => joinpath(get_global(:path_to_data), "crystals"))
    # set the default covalent radii
    set_global(:covalent_radii => DEFAULT_COVALENT_RADII)
    # set the default bonding rules
    set_global(:bonding_rules => bondingrules()) # wanted to precompile, but didn't work
    # set the default atomic masses
    set_global(:atomic_masses => DEFAULT_ATOMIC_MASSES)
    # the set default CPK df_colors
    set_global(:cpk_colors => DEFAULT_CPK_COLORS)
    @debug "Default environment variables" GLOBAL
end


export
    # paths_and_globals.jl
    print_file_paths, set_global, get_global,

    # matter.jl
    Coords, Frac, Cart, Atoms, Charges, wrap!, neutral, net_charge,
    translate_by!, origin,

    # box.jl
    Box, replicate, unit_cube, write_vtk, inside, fractional_coords,
    cartesian_coords,

    # distance.jl
    nearest_image!, distance, overlap, remove_duplicates, pairwise_distances,

    # misc.jl
    read_xyz, write_xyz, read_mol, write_mol2,

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
