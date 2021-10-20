#=
    julia cross_dep.jl PKG

Installs and tests PKG (for the cross-testing CI stage)
=#

import Pkg
ENV["CI_BUILD"] = "true"
Pkg.develop(ARGS[1])
Pkg.test(ARGS[1])
