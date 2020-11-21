module Xtals

using CSV
using DataFrames
using Printf
using LinearAlgebra
using LightGraphs

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

# this runs everytime Xtals is loaded, so if the user changes directory
#   then the path_to_data will change as well
function __init__()
    set_default_file_paths(print_paths=false)
end


include("matter.jl")
include("box.jl")
include("distance.jl")
include("misc.jl")
include("crystal.jl")
include("bonds.jl")
include("repfactors.jl")


export
    # Xtals.jl
    print_file_paths, set_path_to_data,

    # matter.jl
    Coords, Frac, Cart, Atoms, Charges, wrap!, neutral, net_charge,
	translate_by!, origin,

    # box.jl
    Box, replicate, unit_cube, write_vtk, inside, fractional_coords,
	cartesian_coords,

    # distance.jl
    nearest_image!, distance, overlap, remove_duplicates,

    # misc.jl
    read_xyz, read_cpk_colors, write_xyz, read_atomic_masses,

    # crystal.jl
    Crystal, strip_numbers_from_atom_labels!, assign_charges,
    chemical_formula, molecular_weight, crystal_density, write_cif, has_charges,
    apply_symmetry_operations,

    # bonds.jl
    infer_bonds!, write_bond_information, BondingRule, bond_sanity_check,
	remove_bonds!, infer_geometry_based_bonds!, cordero_covalent_atomic_radii

end # module Xtals
