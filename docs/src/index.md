![Xtals.jl](assets/logo.jpg)

A pure-[Julia](https://julialang.org/) package for representation of porous crystals such as metal-organic frameworks (MOFs).

*In development, please contribute, post issues 🐛, and improve!*

## Installation

 1. Download and install the [Julia programming language](https://julialang.org/), v1.3 or higher.

 2. Open the package manager in the Julia REPL (using `]`) and enter:

```
pkg> add Xtals
```

 3. In Julia, load all functions in `Xtals.jl` into the namespace:

```
julia> using Xtals # that's it
```

## Dependencies

Some (non-core) features of `Xtals.jl` require Python libraries.
If these are missing, you may see warnings when the module loads.
For users unfamiliar with configuring Python, the [quick-setup script](https://raw.githubusercontent.com/SimonEnsemble/Xtals.jl/master/quick_setup.jl) will install and build `Xtals.jl` along with a Python environment and the dependencies.
Run it like so:

```
$ julia quick_setup.jl
```

## Tests

Unit tests for the package are available by doing `] test Xtals` in the Julia REPL.
