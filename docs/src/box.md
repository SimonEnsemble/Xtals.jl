```@meta
DocTestSetup = quote
  using Xtals
end
```

# The Spatial Box

Within `Xtals.jl`, the 3D space in which all [`Coords`](@ref) are located is the [`Box`](@ref).
Each [`Crystal`](@ref) has its own [`Box`](@ref), equivalent to the unit cell of a material, containing as attributes the unit cell edge lengths (`a` `b` `c`), crystallographic dihedral angles (`α` `β` `γ`), volume, conversion factors for translating between [`Frac`](@ref)tional and [`Cart`](@ref)esian coordinates, and the reciprocal (Fourier transform) vectors for the Bravais lattice.

## Defining a Box

A [`Box`](@ref) is most conveniently constructed from its basic spatial data (`a` `b` `c` `α` `β` `γ`).
For example, given the unit cell of Co-MOF-74, we can define its [`Box`](@ref):

```julia
a = 26.13173 # Å
b = 26.13173
c = 6.722028
α = π / 2 # radians
β = π / 2
γ = 2 * π / 3
box = Box(a, b, c, α, β, γ)
```

A [`Box`](@ref) may also be defined by providing only the [`Frac`](@ref)tional-to-[`Cart`](@ref)esian conversion matrix:

```julia
box = Box([26.1317 -13.0659 0; 0 22.6307 0; 0 0 6.72203])
```

To quickly get a simple unit-cubic [`Box`](@ref), use the [`unit_cube`](@ref) function.

```julia
@info unit_cube()
#┌ Info: Bravais unit cell of a crystal.
#│       Unit cell angles α = 90.000000 deg. β = 90.000000 deg. γ = 90.000000 deg.
#│       Unit cell dimensions a = 1.000000 Å. b = 1.000000 Å, c = 1.000000 Å
#└       Volume of unit cell: 1.000000 Å³
```

## Transforming Coordinates

Conversions are provided for switching between [`Frac`](@ref)tional and [`Cart`](@ref)esian [`Coords`](@ref) using the [`Box`](@ref) (works for [`Atoms`](@ref) and [`Charges`](@ref), too).

```julia
xtal = Crystal("Co-MOF-74.cif")
Cart(xtal.atoms.coords, xtal.box)
#Cart([-5.496156112249995 7.181391379950001 … 15.131970232450003 2.4686645331000063;
# 22.270234304380295 2.8331425940892103 … 0.7607701110682343 22.13256395706254;
# 1.231811631 0.32198514120000005 … 6.2082409932000004 2.2119953472])
```

## Replicating a Box

For simulations in larger volumes than a single crystallograhic unit cell, the [`Box`](@ref) may be replicated along each or any of the three crystallographic axes.
See [`replicate`](@ref).

```julia
replicated_box = replicate(box, (2, 2, 2))
```

## Exporting a Box

For visualization of the unit cell boundaries, the [`Box`](@ref) may be written out to a `.vtk` file for use in [Visit](https://wci.llnl.gov/simulation/computer-codes/visit/).

```jldoctest; setup=:(box = Box([26.1317 -13.0659 0; 0 22.6307 0; 0 0 6.72203])), output=false
write_vtk(box, "box.vtk"; verbose=true)

# output

See box.vtk
```

# Detailed Docs

```@docs
    Box
    unit_cube
    write_vtk
```
