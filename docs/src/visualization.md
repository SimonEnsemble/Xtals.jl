# Visualizing structures

Crystals can be visualized using [`Bio3DView`](https://github.com/jgreener64/Bio3DView.jl) and either [`Blink`](https://juliagizmos.github.io/Blink.jl/stable/), [`Jupyter`](https://julialang.github.io/IJulia.jl/stable/), or [`Pluto`](https://github.com/fonsp/Pluto.jl).

To include a visualization within a `Pluto` or `Jupyter` notebook, simply call [`view_crystal`](@ref). To launch a visualization window directly from a script or the REPL, use `Blink`, like so:

```julia
using Blink
xtal = Crystal("SBMOF-1.cif")
infer_bonds!(xtal, true)
view_crystal(xtal)
```

## Docs

```@docs
view_crystal
```
