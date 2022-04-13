```@meta
DocTestSetup = quote
  using Xtals
end
```

# Distances

The distance between two [`Atoms`](@ref) in a [`Crystal`](@ref) is central to many operations within `Xtals.jl`.  
The [`distance`](@ref) function calculates the [`Cart`](@ref)esian displacement between the [`Coords`](@ref) ([`Cart`](@ref) or [`Frac`](@ref)) of two points, `i` and `j`, within a given [`Box`](@ref), in units of Å.

```jldoctest distance
xtal = Crystal("IRMOF-1.cif")
distance(xtal.atoms.coords, xtal.box, 1, 2, false)
# output
18.538930020700434
```

The `apply_pbc` argument allows for calculation of distances across the periodic boundaries of the [`Box`](@ref).

```jldoctest distance
distance(xtal.atoms.coords, xtal.box, 1, 2, true)
# output
15.096469110488975
```

[`distance`](@ref) also works on [`Atoms`](@ref) and [`Charges`](@ref).

```jldoctest distance
distance(xtal.atoms, xtal.box, 3, 5, true)
# output
16.90555095103936
```

# Detailed Docs

```@docs
distance
```
