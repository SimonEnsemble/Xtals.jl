```@meta
DocTestSetup = quote
  using Xtals
end
```

# Chemical Bonding

Chemical bonding interactions are represented in the `bonds` attribute of a [`Crystal`](@ref) as a [graph](https://github.com/JuliaGraphs/MetaGraphs.jl).
The nodes of the graph correspond to the [`Crystal`](@ref)'s [`Atoms`](@ref), and the edges of the graph correspond to the bonds (`xtal.bonds`).

## Bonding Rules

`Xtals` uses an array of [`BondingRule`](@ref) structs stored in `rc` for deciding if two [`Atoms`](@ref) are an appropriate distance to be chemically bonded.  
The default rules are based on the [Cordero covalent radii](https://doi.org/10.1039/B801115J), modified based on the work of [Thomas Manz](https://doi.org/10.1039/c9ra07327b).  
Each [`BondingRule`](@ref) is composed of two chemical species symbols and a floating point value, the maximum distance for inferring a bond between the indicated species.

```jldoctest; output=false
BondingRule(:C, :C, 1.77)
# output
BondingRule(:C, :C, 1.77)
```

The global bonding rules may be augmented with [`add_bonding_rules`](@ref) or written to/read from disk with [`write_bonding_rules`](@ref) and [`read_bonding_rules`](@ref).  
The default rules are determined from `rc[:covalent_radii]` at module load, but *are not recalculated upon changes to the covalent radii.*  
If `rc[:covalent_radii]` is altered and new bonding rules should be calculated, the user must do `rc[:bonding_rules] = bondingrules()`.

## Adding Bonds to Crystals

By default, bonding information is not added to a [`Crystal`](@ref). 
Bonds may be inferred at the time of loading structure data by using the `infer_bonds` keyword argument.  
See [`Crystal`](@ref) for more details.

```jldoctest bonds
xtal = Crystal("SBMOF-1.cif", infer_bonds=true, periodic_boundaries=true)
xtal.bonds
# output
{120, 144} undirected Int64 metagraph with Float64 weights defined by :weight (default weight 1.0)
```

The first number is the number of nodes, and the second is the number of edges.  
The `:weight` attribute is not used, and can be ignored.

[`remove_bonds!`](@ref) clears the bonding information from a [`Crystal`](@ref):

```jldoctest bonds
remove_bonds!(xtal)
xtal.bonds
# output
{120, 0} undirected Int64 metagraph with Float64 weights defined by :weight (default weight 1.0)
```

Use [`infer_bonds!`](@ref) to infer plausible bonds using the global bonding rules (or another specified set of rules) in already-loaded [`Crystal`](@ref)s:

```jldoctest bonds
infer_bonds!(xtal, true)
xtal.bonds
# output
{120, 144} undirected Int64 metagraph with Float64 weights defined by :weight (default weight 1.0)
```

## Bonds from Input File

Some chemical information formats, like `.cif` and `.mol`, can store not only the cartesian coordinates of atoms, but also the graph of bonds between atoms in molecules and crystals.  
The `read_bonds_from_file` keyword argument for [`Crystal`](@ref) enables loading these bonds when reading the data.  
[`read_mol`](@ref) also returns bond information.

## Bonds for Atoms

In the case that atomic coordinates are loaded from XYZ format, there will be no unit cell information.  
To infer bonds between atoms in this case, use the [`infer_bonds`](@ref) function:

```jldoctest bonds
atoms = Cart(xtal.atoms, xtal.box) # get atoms in Cartesian coords
bonds = infer_bonds(atoms) # infer bonding graph
# output
{120, 110} undirected Int64 metagraph with Float64 weights defined by :weight (default weight 1.0)
```

## Bond Distances, Vectors, and Angles

Bonds may be labeled with several additional pieces of information.  
The first is the center-to-center distance between the bonded atoms, accessible via [`bond_distance`](@ref):

```jldoctest bonds
bond_distance(xtal, 1, 5)
# output
1.5233240952030063
```

The bond distance is automatically added by [`infer_bonds!`](@ref). 
Applying [`calculate_bond_vectors!`](@ref) (or passing `calculate_vectors=true` to [`infer_bonds!`](@ref)) labels each edge in the bonding graph with a vector, accessible via [`get_bond_vector`](@ref):

```jldoctest bonds
calculate_bond_vectors!(xtal)
get_bond_vector(xtal, 1, 5)
# output
3-element Vector{Float64}:
 -1.2460713575602618
 -0.26441824999999985
  0.8354073616870523
```

While the bond graph itself is undirected, the vectors returned by [`get_bond_vector`](@ref) are directed, so reversing the order of the node indices will flip the vector:

```jldoctest bonds
get_bond_vector(xtal, 5, 1)
# output
3-element Vector{Float64}:
  1.2460713575602618
  0.26441824999999985
 -0.8354073616870523
```

Bond angles are calculated via the dot product: θ = arccos(u⃗•v⃗ / (‖u⃗‖*‖v⃗‖))
To get the angle (in radians) between two bonds, use [`bond_angle`](@ref):

```jldoctest bonds
bond_angle(xtal, 1, 5, 9)
# output
2.089762039489374
```

## Detailed Docs

```@docs
BondingRule
add_bonding_rules
read_bonding_rules
write_bonding_rules
infer_bonds!
infer_bonds
remove_bonds!
write_bond_information
bond_distance
calculate_bond_vectors!
clear_vectors!
get_bond_vector
bond_angle
```
