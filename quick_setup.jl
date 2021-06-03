@info "Running Xtals.jl quick-setup."

import Pkg

@info "Checking Python dependencies..."
try # see if scipy is already available
    using PyCall
    pyimport("scipy")
catch # if not, set up Python environment and get scipy
    @info "Setting up Python dependencies..."
    ENV["PYTHON"]=""
    Pkg.add("PyCall")
    Pkg.build("PyCall")
    Pkg.add("Conda")
    using Conda
    Conda.add("scipy")
end

@info "Checking Xtals..."
try # see if Xtals is already installed
    using Xtals
catch # if not, install it
    @info "Installing Xtals..."
    Pkg.add("Xtals")
    Pkg.build("Xtals")
end

@info "Done!"