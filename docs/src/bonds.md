```@meta
DocTestSetup = quote
  using Xtals
end
```

# Chemical Bonding

Chemical bonding interactions are represented in the `bonds` attribute of a `Crystal` as a [graph](https://github.com/JuliaGraphs/MetaGraphs.jl) where the nodes of the graph correspond to the crystal's `atoms`.

## Bonding Rules

`Xtals` uses an array of [`BondingRule`](@ref) structs stored at [`rc`](@ref) for deciding if two atoms are an appropriate distance to be chemically bonded.  The default rules are based on the [Cordero covalent radii](doi.org/10.1039/B801115J), modified based on the work of [Thomas Manz](doi.org/10.1039/c9ra07327b).  Each [`BondingRule`](@ref) is composed of two chemical species symbols and a floating point value, the maximum distance for inferring a bond between the indicated species.

```jldoctest; output=false
BondingRule(:C, :C, 1.77)
# output
BondingRule(:C, :C, 1.77)
```

The global bonding rules are accessible at [`rc[:bonding_rules]`](@ref) and may be augmented with [`add_bonding_rules`](@ref) or written to/read from disk with [`write_bonding_rules`](@ref) and [`read_bonding_rules`](@ref).  The default rules are determined from `rc[:covalent_radii]` at module load, but *are not recalculated upon changes to the covalent radii.*  If `rc[:covalent_radii]` is altered and new bonding rules should be calculated, the user must do `rc[:bonding_rules] = bondingrules()`.

## Adding Bonds to Crystals

By default, bonding information is not added to a `Crystal`.  Use [`infer_bonds!`](@ref) or [`infer_geometry_based_bonds`](@ref) to automatically infer plausible bonds using the global bonding rules, or another specified set of rules.  [`remove_bonds!`](@ref) clears the bonding information from a `Crystal`.

Bonds may also be inferred at the time of loading crystal data, using the `infer_bonds` keyword argument.  See [`Crystal`](@ref) for more details.

## Bonds from Input File

Some chemical information formats, like `.cif` and `.mol`, can store not only the cartesian coordinates of atoms, but also the graph of bonds between atoms in molecules and crystals.  The `read_bonds_from_file` keyword argument for [`Crystal`](@ref) enables loading these bonds when reading the data.  [`read_mol`](@ref) also returns bond information.

## Detailed Docs

```@docs
BondingRule
add_bonding_rules
read_bonding_rules
write_bonding_rules
infer_bonds!
remove_bonds!
infer_geometry_based_bonds!
write_bond_information
read_mol
```
