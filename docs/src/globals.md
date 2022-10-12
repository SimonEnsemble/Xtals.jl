# Global Variables

In `Xtals.jl`, global variables are stored in a dictionary called `rc`.
The entries in `rc` are used for holding various important pieces of information:

  - `rc[:atomic_masses]` : a dictionary which maps atomic species to their masses in AMU
  - `rc[:cpk_colors]` : a dictionary which maps atomic species to their CPK colors
  - `rc[:bonding_rules]` : an array which lists the maximum bonding distances between each pair of atom types in Å
  - `rc[:covalent_radii]` : a dictionary which maps atomic species to their covalent radii in Å
  - `rc[:paths]` : a dictionary of paths to directories to search for input data

## Paths

The `rc[:paths]` dictionary holds the following:

  - `rc[:paths][:crystals]` : the path to find crystallographic files (`.cif`, `.cssr`) to read in
  - `rc[:paths][:data]` : not directly used by `Xtals.jl`, but important for other packages which use `Xtals`, e.g. [PorousMaterials.jl](https://simonensemble.github.io/PorousMaterials.jl/dev/)

When `Xtals` is first loaded, the paths are set based on the present working directory.
`Xtals` expects to find a folder named `data` in the working directory, and assumes that crystallographic data will be present in `data/crystals`.
To change the location of the crystal inputs, either set `rc[:paths][:crystals]` directly, or use [`set_paths`](@ref):

```julia
rc[:paths][:crystals] = "other_crystals" # crystals are now loaded from "other_crystals"
set_paths("other_data") # crystals are now loaded from "other_data/crystals"
```

# Detailed Docs

```@docs
set_paths
```
