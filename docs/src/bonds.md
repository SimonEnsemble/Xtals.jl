# Chemical Bonding

Chemical bonding interactions are represented in the `bonds` attribute of a `Crystal`
as a graph where the nodes of the graph correspond to the crystal's `atoms`.

```@docs
BondingRule
get_bonding_rules,
set_bonding_rules
add_bonding_rules
read_bonding_rules
write_bonding_rules
infer_bonds!
bond_sanity_check
remove_bonds!
infer_geometry_based_bonds!
write_bond_information
```
