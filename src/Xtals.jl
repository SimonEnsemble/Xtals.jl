module Xtals

using CSV
#using DataFrames
#using Roots # for fzero
#using SpecialFunctions # for erfc
#using StatsBase
#using ProgressMeter
#using Polynomials
#using JLD2
#using Statistics
using Printf
using LinearAlgebra
using LightGraphs
#using Distributed
#using Optim
#using PyCall

# atoms are considered to overlap if this close.
const R²_OVERLAP = 0.1 # Units: Angstrom²

"""
    print_file_paths()

print off paths where Xtals.jl looks for input files and writes output files.
"""
function print_file_paths()
    println("general data folder: ", PATH_TO_DATA)
    println("\tcrystal structures (.cif, .cssr): ", PATH_TO_CRYSTALS)
    println("\tforce field files (.csv): ", PATH_TO_FORCEFIELDS)
    println("\tmolecule input files: ", PATH_TO_MOLECULES)
    println("\tsimulation output files: ", PATH_TO_SIMS)
    println("\tgrids (.cube): ", PATH_TO_GRIDS)
end

"""
    set_path_to_data("../data")

set the `PATH_TO_DATA` variable.
this adjusts the path to crystals, forcefields, molecules, grids, and simulation output files.
"""
function set_path_to_data(ptd::String; print_paths::Bool=true)
    global PATH_TO_DATA = ptd

    global PATH_TO_CRYSTALS = joinpath(PATH_TO_DATA, "crystals")
    global PATH_TO_FORCEFIELDS = joinpath(PATH_TO_DATA, "forcefields")
    global PATH_TO_MOLECULES = joinpath(PATH_TO_DATA, "molecules")
    global PATH_TO_GRIDS = joinpath(PATH_TO_DATA, "grids")
    global PATH_TO_SIMS = joinpath(PATH_TO_DATA, "simulations")

    if print_paths
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
    global PATH_TO_DATA = joinpath(pwd(), "data")

    global PATH_TO_CRYSTALS = joinpath(PATH_TO_DATA, "crystals")
    global PATH_TO_FORCEFIELDS = joinpath(PATH_TO_DATA, "forcefields")
    global PATH_TO_MOLECULES = joinpath(PATH_TO_DATA, "molecules")
    global PATH_TO_GRIDS = joinpath(PATH_TO_DATA, "grids")
    global PATH_TO_SIMS = joinpath(PATH_TO_DATA, "simulations")

    if print_paths
        print_file_paths()
    end
end

# this runs everytime Xtals is loaded, so if the user changes directory
#   then the path_to_data will change as well
function __init__()
    set_default_file_paths(print_paths=false)
end

 # """
 #     set_tutorial_mode()
 #
 # places Xtals in "tutorial mode". it changes the `path_to_data` variable to
 # the directory where the Xtals test data is stored. it can be used to
 # follow examples shown in the readme. it displays a warning so that the user knows
 # they are no longer using their own data.
 # """
 # function set_tutorial_mode()
 #     new_path = joinpath(dirname(pathof(Xtals)), "..", "test", "data")
 #     if ! isdir(new_path)
 #         @error @sprintf("directory for testing data %s does not exist.\nnot entering tutorial mode.\n", new_path)
 #     else
 #         global path_to_data = new_path
 #         global path_to_crystals = joinpath(path_to_data, "crystals")
 #         global path_to_forcefields = joinpath(path_to_data, "forcefields")
 #         global path_to_molecules = joinpath(path_to_data, "molecules")
 #         global path_to_grids = joinpath(path_to_data, "grids")
 #         @warn "Xtals is now in tutorial mode. you have access to the testing data to experiment with Xtals.\nto reset to default file paths, call `set_default_file_paths()`\n"
 #     end
 # end
 #
include("matter.jl")
include("box.jl")
include("distance.jl")
include("misc.jl")
include("crystal.jl")
include("bonds.jl")

"""
	repfactors = replication_factors(unitcell, cutoffradius)
Find the replication factors needed to make a supercell big enough to fit a sphere with the specified cutoff radius.
In Xtals.jl, rather than replicating the atoms in the home unit cell to build the supercell that
serves as a simulation box, we replicate the home unit cell to form the supercell (simulation box) in a for loop.
This function ensures enough replication factors such that the nearest image convention can be applied.
A non-replicated supercell has 1 as the replication factor in each dimension (`repfactors = (1, 1, 1)`).
# Arguments
- `unitcell::Box`: The unit cell of the framework
- `cutoff_radius::Float64`: Cutoff radius beyond which we define the potential energy to be zero (units: Angstrom)
# Returns
- `repfactors::Tuple{Int, Int, Int}`: The replication factors in the a, b, c directions
"""
function replication_factors(unitcell::Box, cutoff_radius::Float64)
	# Unit vectors used to transform from fractional coordinates to cartesian coordinates. We'll be
	a = unitcell.f_to_c[:, 1]
	b = unitcell.f_to_c[:, 2]
	c = unitcell.f_to_c[:, 3]

	n_ab = cross(a, b)
	n_ac = cross(a, c)
	n_bc = cross(b, c)

	# c0 defines a center in the unit cell
	c0 = [a b c] * [.5, .5, .5]

	rep = [1, 1, 1]

	# Repeat for `a`
	# |n_bc ⋅ c0|/|n_bc| defines the distance from the end of the supercell and the center. As long as that distance is less than the cutoff radius, we need to increase it
	while abs(dot(n_bc, c0)) / norm(n_bc) < cutoff_radius
		rep[1] += 1
		a += unitcell.f_to_c[:,1]
		c0 = [a b c] * [.5, .5, .5]
	end

	# Repeat for `b`
	while abs(dot(n_ac, c0)) / norm(n_ac) < cutoff_radius
		rep[2] += 1
		b += unitcell.f_to_c[:,2]
		c0 = [a b c] * [.5, .5, .5]
	end

	# Repeat for `c`
	while abs(dot(n_ab, c0)) / norm(n_ab) < cutoff_radius
		rep[3] += 1
		c += unitcell.f_to_c[:,3]
		c0 = [a b c] * [.5, .5, .5]
	end

	return (rep[1], rep[2], rep[3])::Tuple{Int, Int, Int}
end

export
    # Xtals.jl
    print_file_paths, set_path_to_data,

    # matter.jl
    Coords, Frac, Cart, Atoms, Charges, wrap!, neutral, net_charge, translate_by!, origin,

    # box.jl
    Box, replicate, unit_cube, write_vtk, inside, fractional_coords, cartesian_coords,

    # distance.jl
    nearest_image!, distance, overlap, remove_duplicates,

    # misc.jl
    read_xyz, read_cpk_colors, write_xyz, read_atomic_masses,

    # crystal.jl
    Crystal, strip_numbers_from_atom_labels!, assign_charges,
    chemical_formula, molecular_weight, crystal_density, write_cif, has_charges,
    apply_symmetry_operations,

    # bonds.jl
    infer_bonds!, write_bond_information, BondingRule, bond_sanity_check, remove_bonds!,
    infer_geometry_based_bonds!, cordero_covalent_atomic_radii
end
