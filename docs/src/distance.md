```@meta
DocTestSetup = quote
  using Xtals
end
```

# Distances

The distance between two `Atoms` in a `Crystal` is central to many operations
within `Xtals.jl`.  The `distance` function calculates the `Cart`esian
displacement between the `Coords` (`Cart` or `Frac`) of two points, `i` and `j`,
within a given `box`, in units of Å.

```jldoctest distance
xtal = Crystal("IRMOF-1.cif")
distance(xtal.atoms.coords, xtal.box, 1, 2, false)
# output
18.538930020700434
```

The `apply_pbc` argument allows for calculation of distances
across the periodic boundaries of the `box`.

```jldoctest distance
distance(xtal.atoms.coords, xtal.box, 1, 2, true)
# output
15.096469110488975
```

`distance` also works on `Atoms` and `Charges`.

```jldoctest distance
distance(xtal.atoms, xtal.box, 3, 5, true)
# output
16.90555095103936
```

# docs

```@docs
    distance
```
