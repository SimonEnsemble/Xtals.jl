# XtalsPyTools.jl

`Xtals` is a pure-Julia package, but we have also written some high-level tools with it that depend on Python packages.
These additional methods are found in [`XtalsPyTools.jl`](https://github.com/SimonEnsemble/XtalsPyTools.jl).

Long-term, we would like refactor these tools to depend only on Julia packages, and fold them back in to `Xtals.jl`.
Please [open an issue](https://github.com/SimonEnsemble/Xtals.jl/issues) or make a pull request if you know how this can be done at present!

## Installation

To gain access to these tools, add the package in the REPL package manager:

`pkg> add XtalsPyTools`

The methods in this package require `scipy` and/or `pymatgen` to be installed.
To configure these dependencies automatically in a Julia-managed environment, use the [`quick_setup.jl` script](https://github.com/SimonEnsemble/XtalsPyTools.jl/blob/master/quick_setup.jl):

`julia quick_setup.jl`

## Advanced Bond Inference

In addition to the [`infer_bonds!`](@ref) function, there is the [`infer_geometry_based_bonds!`](@ref) function, which uses `scipy` to compute the Voronoi tesselation of the atoms in a crystal.
Bonds will then only be inferred between atoms which occupy neighboring Voronoi cells.

## Finding the Primitive Unit Cell

Many crystal structures found in databases are not the minimal representation, a.k.a. primitive unit cell.
The primitive cell can be obtained easily by calling [`primitive_cell`](@ref).
This function requires the `pymatgen` python dependency.

```@docs
infer_geometry_based_bonds!
primitive_cell
```
