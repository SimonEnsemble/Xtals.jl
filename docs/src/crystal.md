```@meta
DocTestSetup = quote
  using Xtals
end
```

# Crystals

`Xtals.jl` maintains a data structure `Crystal` that stores information about a crystal structure file.

## Reading in a crystal structure file

Currently, the crystal structure file reader accepts `.cif` and `.cssr` file formats. `Xtals.jl` looks for the crystal structure files in `rc[:paths][:crystals]` which is by default `./data/crystals` (relative to `pwd()` at module loading). By typing `rc[:paths][:crystals] = "my_crystal_dir"`, `Xtals.jl` now looks for the crystal structure file in `my_crystal_dir`.
The files can be read as:

```jldoctest crystal; output=false
xtal = Crystal("IRMOF-1.cif")       # The crystal reader stores the information in xtal
xtal.name                           # The name of the crystal structure file
xtal.box                            # The unit cell information
xtal.atoms                          # The atom coordinates (in fractional space) and the atom identities
xtal.charges                        # The charge magnitude and coordinates (in fractional space)
xtal.bonds                          # Bonding information in the structure. By default this is an empty graph,
                                    #  but use `read_bonds_from_file=true` argument in `Crystal` to read from crystal structure file
xtal.symmetry                       # Symmetry information of the crystal. By default converts the symmetry to P1 symmetry.
                                    #  Use `convert_to_p1=false` argument in `Crystal` to keep original symmetry
# output
Xtals.SymmetryInfo(["x"; "y"; "z";;], "P1", true)
```

## Fixing atom species labels

Often, the atoms species are appended by numbers. This messes with the internal workings of `Xtals.jl`.
To circumvent this problem, the function `strip_numbers_from_atom_labels!(xtal)` removes the appended numbers.
It is important to use this function prior to GCMC or Henry coefficient calculations and bond inference operations.

```jldoctest crystal
xtal = Crystal("IRMOF-1.cif", species_col=["_atom_site_label"])
xtal.atoms.species[1]
# output
:Zn1
```
```jldoctest crystal
strip_numbers_from_atom_labels!(xtal)
xtal.atoms.species[1]
# output
:Zn
```

## Converting the coordinates to cartesian space

The coordinates of the crystals are stored in fractional coordinates. If one needs to analyze the cartesian coordinates of the crystal,
that can be done by using the unit cell information of the crystal.
```julia
xtal.atoms.coords.xf                                    # array of fractional coordinates
cart_coords = xtal.box.f_to_c * xtal.atoms.coords.xf    # array of cartesian coordinates
```

## Creating a super cell

For many simulations, one needs to replicate the unit cell multiple times to create a bigger super cell.

```jldoctest crystal
super_xtal = replicate(xtal, (2,2,2))       # Replicates the original unit cell once in each dimension
xtal.atoms.n, super_xtal.atoms.n
# output
(424, 3392)
```

## Finding other properties

```julia
rho = crystal_density(xtal)         # Crystal density of the crystal in kg/m^2
mw = molecular_weight(xtal)         # The molecular weight of the unit cell in amu
formula = empirical_formula(xtal)    # The irreducible chemical formula of the crystal
```

## Assigning new charges

If the crystal structure file does not contains partial charges, we provide methods to assign new charges to the crystal

```jldoctest crystal; output=false
species_to_charge = Dict(:Zn => 2.0, :C => 0.0, :H => 0.0, :O => -0.61538)  # This method assigns a static charge to atom species
charged_xtal = assign_charges(xtal, species_to_charge, 1e-3)                # This function creates a new charged `Crystal` object.
                                                                            #   The function checks for charge neutrality with a tolerance of 1e-3
new_charges = Charges(zeros(xtal.atoms.n), xtal.atoms.coords)
other_charged_xtal = Crystal(xtal.name, xtal.box, xtal.atoms,               # Here we create a new `Charges` object using an array of new charges.
                             new_charges, xtal.bonds, xtal.symmetry)        #   The number of charges in the array has to be equal to the number of atoms
                                                                            #   and finally a new `Crystal` object is manually created
# output
Name: IRMOF-1.cif
Bravais unit cell of a crystal.
	Unit cell angles α = 90.000000 deg. β = 90.000000 deg. γ = 90.000000 deg.
	Unit cell dimensions a = 25.832000 Å. b = 25.832000 Å, c = 25.832000 Å
	Volume of unit cell: 17237.492730 Å³

	# atoms = 424
	# charges = 424
	chemical formula: C₁₉₂H₉₆O₁₀₄Zn₃₂
	space Group: P1
	symmetry Operations:
		'x, y, z'
	bonding graph:
		# vertices = 424
		# edges = 0
```

## Writing crystal files

We provide methods to write both `.xyz` and `.cif` files

```jldoctest crystal; output=false
write_cif(xtal, "my_new_cif_file.cif")      # Stored in the current directory
write_xyz(xtal, "my_new_xyz_file.xyz")      # stored in the current directory
# output

```


# Detailed docs

```@docs
    Crystal
    SymmetryInfo
    strip_numbers_from_atom_labels!
    replicate
    molecular_weight
    crystal_density
    empirical_formula
    assign_charges
    write_cif
    write_xyz
    read_xyz
    read_mol
```
