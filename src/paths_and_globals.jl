"""
    global_dict()

Returns a copy of the global variable dictionary for the Xtals.jl ecosystem.
Used only for reference; to change values, use `@ref(set_global)`
"""
function global_dict()
    return copy(GLOBAL)
end


"""
    get_global(key::Symbol)

Returns the value for a global variable in the Xtals.jl ecosystem.

# Arguments
- `key::Symbol` : the name of the variable to retrieve
"""
function get_global(key::Symbol)
    return GLOBAL[key]
end


"""

Sets the value for a global variable in the Xtals.jl ecosystem.
If the key is already in use, the value will be overwritten.  If not, it will be added.

# Arguments
- `key_val_pair::Pair{Symbol,Any}` : a key => value pair defining the value to set
"""
function set_global(key_val_pair::Pair)
    merge!(GLOBAL, Dict(key_val_pair))
end


"""
    print_file_paths()

print off paths where Xtals.jl looks for input files and writes output files.
"""
function print_file_paths()
    println("general data folder: ", get_global(:path_to_data))
    println("crystal structures (.cif, .cssr): ", get_global(:path_to_crystals))
end
