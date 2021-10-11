#=
    julia cross_dep.jl PKG

Installs and tests PKG (for the cross-testing CI stage)
=#

import Pkg
Pkg.dev(ARGS[1])
Pkg.test(ARGS[1])
