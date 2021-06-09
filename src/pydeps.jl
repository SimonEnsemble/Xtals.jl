# Python dependencies.
# Pairs are rc_key => library
PYDEPS = Dict(
    :scipy => "scipy.spatial",
    :pymatgen => "pymatgen.io.cif"
)

# Attempts to load a python dependency.
# Calling functions which require these dependencies without successfully loading them will throw exceptions.
function load_pydep(pydep::String)
    try
        return pyimport(pydep)
    catch
        return nothing
    end
end


# Checks for ability to load a python dependency. Intended for precompilation warnings.
function check_pydep(pydep::Pair)
    if isnothing(load_pydep(pydep[2]))
        @warn "Error loading $(pydep[1]). Some functionaltiy may be missing."
        return missing
    end
end


# Called at module initialization to spin up deps
function init_pydeps()
    for pydep ∈ PYDEPS
        rc[pydep[1]] = load_pydep(pydep[2])
    end
end


# Pre-compilation checks
for pydep ∈ PYDEPS
    check_pydep(pydep)
end
