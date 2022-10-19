# global cpk color dictionary
# atom_symbol => (R, G, B)
# built from PeriodicTable cpk_hex entires
DEFAULT_CPK_COLORS = Dict(Symbol(elements[i].symbol) => ([parse(Int, elements[i].cpk_hex[x], base=16) for x in [[2,3],[4,5],[6,7]]]...,) for i in eachindex(elements) if length(elements[i].cpk_hex) == 7)
