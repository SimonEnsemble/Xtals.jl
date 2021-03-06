# Chemical Bonding

Chemical bonding interactions are represented in the `bonds` attribute of a `Crystal`
as a [graph](https://github.com/JuliaGraphs/MetaGraphs.jl) where the nodes of the
graph correspond to the crystal's `atoms`.

## Bonding Rules

`Xtals` uses a global array of [`BondingRule`](@ref) structs for deciding if two atoms
are an appropriate distance to be chemically bonded.  The default rules are based
on the Cordero covalent radii.  Each [`BondingRule`](@ref) is composed of two chemical
species symbols and two floating point values, the minimum and maximum distances for
inferring a bond between the indicated species.

```julia
BondingRule(:C, :C, 1.02, 1.77)
```

The global bonding rules are available via [`get_bonding_rules`](@ref) and may be
appended with [`add_bonding_rules`](@ref), overwritten with [`set_bonding_rules`](@ref),
or written to/read from disk with [`write_bonding_rules`](@ref) and
[`read_bonding_rules`](@ref).

## Adding Bonds to Crystals

By default, bonding information is not added to a `Crystal`.  Use [`infer_bonds!`](@ref)
or [`infer_geometry_based_bonds`](@ref) to automatically infer plausible bonds using
the global bonding rules, or another specified set of rules.  [`remove_bonds!`](@ref)
clears the bonding information from a `Crystal`.

Bonds may also be inferred at the time of loading crystal data, using the `infer_bonds`
keyword argument.  See [`Crystal`](@ref) for more details.

## Bonds from Input File

Some chemical information formats, like `.cif` and `.mol`, can store not only the
cartesian coordinates of atoms, but also the graph of bonds between atoms in molecules
and crystals.  The `read_bonds_from_file` keyword argument for [`Crystal`](@ref)
enables loading these bonds when reading the data.  [`read_mol`](@ref) also returns
bond information.

## Detailed Docs

```@docs
BondingRule
Xtals.bondingrules
Xtals.get_covalent_radii
get_bonding_rules
set_bonding_rules
add_bonding_rules
read_bonding_rules
write_bonding_rules
infer_bonds!
remove_bonds!
infer_geometry_based_bonds!
write_bond_information
read_mol
```
