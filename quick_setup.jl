@info "Running Xtals.jl quick-setup."

import Pkg
Pkg.build()

@info "Setting up Python dependencies..."
ENV["PYTHON"]=""
Pkg.add("PyCall")
Pkg.build("PyCall")
Pkg.add("Conda")
using Conda
Conda.add("scipy")

@info "Setting up Xtals..."
Pkg.add("Xtals")
Pkg.build("Xtals")

@info "Done!"