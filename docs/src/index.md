# Xtals.jl

A pure-[Julia](https://julialang.org/) package for representation of porous
crystals such as metal-organic frameworks (MOFs).

*In development, please contribute, post issues 🐛, and improve!*

## Installation

1. Download and install the [Julia programming language](https://julialang.org/),
 v1.6 or higher.

2. If the system does not have Python and scipy installed, run the following:

```
julia> import Pkg
julia> Pkg.build()
julia> ENV["PYTHON"]=""
julia> Pkg.build("PyCall")
julia> Pkg.add("Conda")
julia> using Conda
julia> Conda.add("scipy")
```

3. In Julia, open the package manager (using `]`) and enter the following:

```
pkg> add Xtals
```

4. In Julia, load all functions in `Xtals.jl` into the namespace:

```
julia> using Xtals # that's it
```

## Tests
Run the tests in the script `tests/runtests.jl` manually or by `] test Xtals` in the Julia REPL.
