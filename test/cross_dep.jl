#=
    julia cross_dep.jl PKG

Installs and tests PKG (for the cross-testing CI stage)
=#

import Pkg
ENV["CI_BUILD"] = "true" # PoreMatMod needs this; may be useful elsewhere, too
Pkg.develop(ARGS[1])
Pkg.test(ARGS[1])
