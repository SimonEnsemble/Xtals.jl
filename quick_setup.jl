import Pkg


# tries to load a Python dependency; on failure, adds the dependency via Conda
function check_add_dep(pkg; channel="")
    try
        @info "Checking dependency: $pkg"
        pyimport(pkg)
    catch
        @info "Installing $pkg..."
        Conda.add(pkg, channel=channel)
        pyimport(pkg)
        @info "$pkg install verified."
    end
end


@info "Setting up Python environment..."
ENV["PYTHON"]=""
Pkg.add("PyCall")
Pkg.build("PyCall")
using PyCall
pyimport("sys")
@info "PyCall verified."

# check for Conda; if not found, install it
try
    @info "Checking Conda."
    using Conda
    @info "Conda verified."
catch
    @info "Installing Conda..."
    Pkg.add("Conda")
    using Conda
    @info "Conda verified."
end

# check for Xtals; if not found, add it
try
    @info "Checking Xtals."
    using Xtals
    @info "Xtals verified."
catch # if not, install it
    @info "Installing Xtals..."
    Pkg.add("Xtals")
    Pkg.build("Xtals")
    using Xtals
    @info "Xtals install verified."
end

# check for XtalsPyTools; if not found, add it
try
    @info "Checking XtalsPyTools."
    using XtalsPyTools
    @info "XtalsPyTools verified."
catch # if not, install it
    @info "Installing XtalsPyTools..."
    Pkg.add("XtalsPyTools")
    Pkg.build("XtalsPyTools")
    using XtalsPyTools
    @info "XtalsPyTools install verified."
end

# check the deps, add if missing
check_add_dep("scipy")
check_add_dep("pymatgen", channel="conda-forge")

@info "Setup complete!"
