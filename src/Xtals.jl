module Xtals

using CSV, DataFrames, Printf, LinearAlgebra, LightGraphs, PyCall, MetaGraphs

# atoms are considered to overlap if this close.
const R²_OVERLAP = 0.1 # Units: Angstrom²

"""
    print_file_paths()

print off paths where Xtals.jl looks for input files and writes output files.
"""
function print_file_paths()
    println("general data folder: ", PATH_TO_DATA)
    println("crystal structures (.cif, .cssr): ", PATH_TO_CRYSTALS)
end

"""
    set_path_to_data("../data/")

set the `PATH_TO_DATA` or `PATH_TO_CRYSTALS` variable. by default, sets
`PATH_TO_DATA` to `path`, and `PATH_TO_CRYSTALS` to `path/crystals`.

# Arguments
- `path::String` The path to use for setting the environment variable.
- `relpath_xtals::Bool` Specify `true` to update path to crystals relative to new data path.
- `print::Bool` Specify `true` to print path variables.
"""
function set_path_to_data(path::String; relpath_xtals::Bool=false, print::Bool=false)
    global PATH_TO_DATA = path
    if relpath_xtals
        set_path_to_crystals(joinpath(PATH_TO_DATA, "crystals"))
    end
    if print
        print_file_paths()
    end
end

"""
    set_path_to_crystals("../other_crystals/")

set `Xtals.PATH_TO_CRYSTALS`.

# Arguments
- `path::String` The path to use for setting the environment variable.
- `print::Bool` Specify `true` to print path variables.
"""
function set_path_to_crystals(path::String; print::Bool=false)
    global PATH_TO_CRYSTALS = path
    if print
        print_file_paths()
    end
end

"""
    set_default_file_paths(print_paths=true)

sets the default paths for where input files and some output files are stored.
to see current set up, call [`print_file_paths`](@ref)
"""
function set_default_file_paths(;print_paths::Bool=true)
    # this is the main directory where crystal structures, forcefields, and molecules data is stored
    set_path_to_data(joinpath(pwd(), "data"), relpath_xtals=true)
    if print_paths
        print_file_paths()
    end
end


include("matter.jl")
include("box.jl")
include("crystal.jl")
include("distance.jl")
include("misc.jl")
include("cordero.jl")
include("bonds.jl")
include("repfactors.jl")
include("atomic_masses.jl")
include("cpk_colors.jl")


# runs every time the module is imported
function __init__()
    # if the user changes directory, path_to_data will change as well
    set_default_file_paths(print_paths=false)
    # set the default bonding rules from internal cordero parameters
    set_bonding_rules(bondingrules()) # wanted to precompile, but didn't work
end


export
    # Xtals.jl
    print_file_paths, set_path_to_data, set_path_to_crystals,

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

    # cpk_colors.jl
    read_cpk_colors, get_cpk_colors,

    # atomic_masses.jl
    read_atomic_masses, get_atomic_masses,

    # crystal.jl
    Crystal, strip_numbers_from_atom_labels!, assign_charges,
    chemical_formula, molecular_weight, crystal_density, write_cif, has_charges,
    apply_symmetry_operations,

    # bonds.jl
    infer_bonds!, write_bond_information, BondingRule, bond_sanity_check,
    remove_bonds!, infer_geometry_based_bonds!, get_bonding_rules,
    set_bonding_rules, read_bonding_rules, write_bonding_rules, add_bonding_rules,
    drop_cross_pb_bonds!

end # module Xtals
