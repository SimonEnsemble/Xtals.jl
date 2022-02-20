"""
    bonding_rule = BondingRule(:Ca, :O, 2.0)
    bonding_rules = [BondingRule(:H, :*, 1.2),
                     BondingRule(:*, :*, 1.9)]

A rule for determining if two atoms within a crystal are bonded.

# Attributes
- `species_i::Symbol`: One of the atoms types for this bond rule
- `species_j::Symbol`: The other atom type for this bond rule
- `max_dist`: The maximum distance between the atoms for bonding to occur
"""
struct BondingRule
    species_i::Symbol
    species_j::Symbol
    max_dist::Float64
end

BondingRule(species_i::String, species_j::String, max_dist::String) = BondingRule(
    Symbol.([species_i, species_j])..., parse(Float64, max_dist))


"""
    bond_rules = bondingrules()
    bond_rules = bondingrules(covalent_radii=covalent_radii(), pad=0.05)

Returns a set of bonding rules based on the given Cordero parameters and tolerances.

# Arguments

- `cov_rad::Union{Dict{Symbol, Float64}, Nothing}`: Covalent radii. See [`covalent_radii()`](@ref)
- `pad::Float`: Amount to pad covalent radii for bonding interactions.

# Returns
- `bondingrules::Array{BondingRule, 1}`: Bonding rules from the specified covalent radii.`
"""
function bondingrules(;covalent_radii::Dict{Symbol,Float64}=rc[:covalent_radii], pad::Float64=0.)::Array{BondingRule}
    bond_rules = BondingRule[]
    # loop over parameterized atoms
    for (i, atom_i) in enumerate(keys(covalent_radii))
        # make rules for the atom with every other atom (and itself)
        for (j, atom_j) in enumerate(keys(covalent_radii))
            if j < i
                continue # already did this atom in outer loop (don't duplicate)
            end
            radii_sum = covalent_radii[atom_i] + covalent_radii[atom_j]
            max_dist = radii_sum + pad
            push!(bond_rules, BondingRule(atom_i, atom_j, max_dist))
        end
    end
    return bond_rules
end


"""
    write_bonding_rules("file.csv")
Writes bonding rules to a CSV file that can be loaded with [`read_bonding_rules`](@ref)
# Arguments
- `filename::String` : The name of the output file
- `bonding_rules::Array{BondingRule}` : (Optional) The rules to write to file. If not specified, the global rules are written.
"""
function write_bonding_rules(filename::String, bonding_rules::Array{BondingRule}=rc[:bonding_rules])
    f = open(filename, "w")
    for r ∈ bonding_rules
        @printf(f, "%s,%s,%f\n", r.species_i, r.species_j, r.max_dist)
    end
    close(f)
end


"""
    read_bonding_rules("file.csv")
Reads a CSV file of bonding rules and returns a BondingRule array.
# Arguments
- `filename::String` : name of file in data directory containing bonding rules
# Returns
`rules::Array{BondingRule}` : the bonding rules read from file
"""
function read_bonding_rules(filename::String)::Array{BondingRule}
    rules = BondingRule[]
    open(filename) do input_file
        for line in eachline(input_file)
            push!(rules, BondingRule(String.(split(line, ","))...))
        end
    end
    return rules
end


"""
    add_bonding_rules(bonding_rules)
Adds `bonding_rules` to the beginning of the global bonding rules list
# Arguments
- `bonding_rules::Array{BondingRule}` : the array of bonding rules to add
"""
function add_bonding_rules(bonding_rules::Array{BondingRule})
    br = rc[:bonding_rules]
    rc[:bonding_rules] = vcat(bonding_rules, br)
end
add_bonding_rules(bonding_rule::BondingRule) = add_bonding_rules([bonding_rule])


# for pretty-printing the bonding rules
function show(io::IO, bonding_rules::Array{BondingRule})
    for r in bonding_rules
        println("%s\t%s\t%.3f\t%.3f", r.species_i, r.species_j, r.max_dist)
    end
end

"""
    bonded, bond_length = is_bonded(atoms, i, j, bonding_rules)

Checks to see if atoms `i` and `j` in `crystal` are bonded according to the `bonding_rules`.

# Arguments
- `atoms::Atoms{Cart}`: The atoms that may be bonded
- `i::Int`: Index of the first atom
- `j::Int`: Index of the second atom
- `bonding_rules::Array{BondingRule, 1}`: The array of bonding rules that will
    be used to fill the bonding information. They are applied in the order that
    they appear.

# Returns
- `bonded_flag::Bool`: true iff atoms `i` and `j` are bonded according to `bonding_rules`
- `r::Float64`: the distance between atomic centers.
"""
function is_bonded(atoms::Atoms{Cart}, i::Int64, j::Int64, bonding_rules::Array{BondingRule, 1})
    species_i = atoms.species[i]
    species_j = atoms.species[j]
    # compute vector between the two atoms.
    dxc = atoms.coords.x[:, i] - atoms.coords.x[:, j]
    # compute Euclidean distance
    r = norm(dxc)
    
    bonded_flag = loop_over_bonding_rules(bonding_rules, species_i, species_j, r)

    return bonded_flag, r
end

"""
    bonded = loop_over_bonding_rules(bonding_rules, species_i, species_j, r)
"""
function loop_over_bonding_rules(bonding_rules::Vector{BondingRule}, species_i::Symbol, species_j::Symbol, r::Float64)
    # loop over possible bonding rules
    for br in bonding_rules
        # determine if the atom species correspond to the species in `bonding_rules`
        species_match = false
        if br.species_i == :* && br.species_j == :*
            species_match = true
        elseif br.species_i == :* && (species_i == br.species_j || species_j == br.species_j)
            species_match = true
        elseif br.species_j == :* && (species_i == br.species_i || species_j == br.species_j)
            species_match = true
        elseif (species_i == br.species_i && species_j == br.species_j) || (species_j == br.species_i && species_i == br.species_j)
            species_match = true
        end
        if species_match
            # determine if the atoms are close enough to bond
            if br.max_dist > r
                return true
            else
                # return (don't continue!) b/c we found relevant bonding rule, don't apply others
                return false
            end
        end
    end
    # no bonding rule was applied; therefore not bonded.
    @warn "No bonding rule for $(species_i)-$(species_j) interaction."
    return false
end


"""
    are_atoms_bonded = is_bonded(crystal, i, j, bonding_rules=[BondingRule(:H, :*, 1.2), BondingRule(:*, :*, 1.9)],
                                 include_bonds_across_periodic_boundaries=true)

Checks to see if atoms `i` and `j` in `crystal` are bonded according to the `bonding_rules`.

# Arguments
- `crystal::Crystal`: The crystal that bonds will be added to
- `i::Int`: Index of the first atom
- `j::Int`: Index of the second atom
- `bonding_rules::Array{BondingRule, 1}`: The array of bonding rules that will
    be used to fill the bonding information. They are applied in the order that
    they appear.
- `include_bonds_across_periodic_boundaries::Bool`: Whether to check across the
    periodic boundary when calculating bonds

# Returns
- `bonded_flag::Bool`: true iff atoms `i` and `j` are bonded according to `bonding_rules`
- `r::Float64`: the periodic distance between atomic centers. For non-bonded pairs, `r=NaN`.
- `cross_boundary_flag::Bool`: true iff the bond is across a periodic boundary.
"""
function is_bonded(crystal::Crystal, i::Int64, j::Int64, bonding_rules::Array{BondingRule, 1};
                   include_bonds_across_periodic_boundaries::Bool=true)
    species_i = crystal.atoms.species[i]
    species_j = crystal.atoms.species[j]
    # compute vector between the two atoms.
    dxf = crystal.atoms.coords.xf[:, i] - crystal.atoms.coords.xf[:, j]
    # apply nearest image convention
    #  on the way, determine if the edge is in the home unit cell or across periodic boundary.
    cross_boundary_flag = false
    if include_bonds_across_periodic_boundaries
        for i in eachindex(dxf)
            @inbounds if abs(dxf[i]) > 0.5
                cross_boundary_flag = true
                @inbounds dxf[i] -= sign(dxf[i])
            end
        end
    end
    # convert dxf to Cartesian coordinates; compute periodic distance
    r = norm(crystal.box.f_to_c * dxf)
    
    bonded = loop_over_bonding_rules(bonding_rules, species_i, species_j, r)

    r = bonded ? r : NaN

    return bonded, r, cross_boundary_flag
end


"""
    remove_bonds!(crystal)

Remove all bonds from a crystal structure, `crystal::Crystal`.
"""
function remove_bonds!(crystal::Crystal)
    while ne(crystal.bonds) > 0
        rem_edge!(crystal.bonds, collect(edges(crystal.bonds))[1].src,
            collect(edges(crystal.bonds))[1].dst)
    end
    clear_vectors!(crystal)
end


"""
    bonds = infer_bonds(atoms; bonding_rules=rc[:bonding_rules])

Returns a `MetaGraph` encoding the chemical bonds between atoms.

# Arguments
- `atoms::Atoms{Cart}`: the atoms to bond
- `bonding_rules::Vector{BondingRule}`: the bonding rules to use
"""
function infer_bonds(atoms::Atoms{Cart}; bonding_rules=rc[:bonding_rules])::MetaGraph
    bonds = MetaGraph(atoms.n)

    for i in 1:atoms.n
        # loop over every unique pair of atoms
        for j in i+1:atoms.n
            bonding_flag, r = is_bonded(atoms, i, j, bonding_rules)
            if bonding_flag
                add_edge!(bonds, i, j)
                set_prop!(bonds, i, j, :distance, r)
                set_prop!(bonds, i, j, :cross_boundary, false)
            end
        end
    end

    return bonds
end


"""
    infer_bonds!(crystal, include_bonds_across_periodic_boundaries; bonding_rules=rc[:bonding_rules])

Populate the bonds in the crystal object based on the bonding rules. If a pair doesn't have a suitable rule then they will not be considered bonded.

`:*` is considered a wildcard and can be substituted for any species. It is a good idea to include a bonding rule between two `:*` to allow any atoms to bond as long as they are close enough.

The bonding rules are hierarchical, i.e. the first bonding rule takes precedence over the latter ones.

# Arguments
- `crystal::Crystal`: The crystal that bonds will be added to.
- `include_bonds_across_periodic_boundaries::Bool`: Whether to check across the periodic boundary when calculating bonds.
- `bonding_rules::Array{BondingRule, 1}`: The array of bonding rules that will be used to fill the bonding information. They are applied in the order that they appear. `rc[:bonding_rules]` will be used if none provided.
- `calculate_vectors::Bool`: Optional. Set `true` to annotate all edges in the `bonds` graph with vector information.
"""
function infer_bonds!(crystal::Crystal, include_bonds_across_periodic_boundaries::Bool;
                      bonding_rules::Array{BondingRule, 1}=rc[:bonding_rules], calculate_vectors::Bool=false)
    @assert ne(crystal.bonds) == 0 @sprintf("The crystal %s already has bonds. Remove them with the `remove_bonds!` function before inferring new ones.", crystal.name)
    # loop over every atom
    for i in 1:crystal.atoms.n
        # loop over every unique pair of atoms
        for j in i+1:crystal.atoms.n
            bonding_flag, r, cross_boundary_flag = is_bonded(crystal, i, j, bonding_rules;
                                                             include_bonds_across_periodic_boundaries=include_bonds_across_periodic_boundaries)
            if bonding_flag
                add_edge!(crystal.bonds, i, j)
                set_prop!(crystal.bonds, i, j, :distance, r)
                set_prop!(crystal.bonds, i, j, :cross_boundary, cross_boundary_flag)
            end
        end
    end
    if calculate_vectors
        calculate_bond_vectors!(crystal)
    end
    bond_sanity_check(crystal)
end


"""
    sane_bonds = bond_sanity_check(crystal)

Run sanity checks on `crystal.bonds`.
* is the bond graph fully connected? i.e. does every vertex (=atom) in the bond graph have at least one edge?
* each hydrogen can have only one bond
* each carbon can have a maximum of four bonds

if sanity checks fail, refer to [`write_bond_information`](@ref) to write a .vtk to visualize the bonds.

Print warnings when sanity checks fail.
Return `true` if sanity checks pass, `false` otherwise.
"""
function bond_sanity_check(crystal::Crystal)::Bool
    for a = 1:crystal.atoms.n
        ns = neighbors(crystal.bonds, a)
        # is the graph fully connected?
        if length(ns) == 0
            @warn "atom $a = $(crystal.atoms.species[a]) in $(crystal.name) is not bonded to any other atom."
            return false
        end
        # does hydrogen have only one bond?
        if (crystal.atoms.species[a] == :H) && (length(ns) > 1)
            @warn "hydrogen atom $a in $(crystal.name) is bonded to more than one atom!"
            return false
        end
        # does carbon have greater than four bonds?
        if (crystal.atoms.species[a] == :C) && (length(ns) > 4)
            @warn "carbon atom $a in $(crystal.name) is bonded to more than four atoms!"
            return false
        end
    end
    return true
end


"""
    write_bond_information(crystal, filename)
    write_bond_information(crystal, center_at_origin=false)
    write_bond_information(xtal, filename, :cross_boundary => p -> p, "bonds.vtk") # cross-boundary bonds only
    write_bond_information(xtal, filename, :distance => d -> d < 1.0, "bonds.vtk") # distance less than 1.0

Writes the bond information from a crystal to the selected filename.

# Arguments
- `crystal::Crystal`: The crystal to have its bonds written to a vtk file
- `filename::String`: The filename the bond information will be saved to. If left out, will default to crystal name.
- `center_at_origin::Bool`: (optional) center the coordinates at the origin of the crystal
- `bond_filter::Pair{Symbol, Function}`: (optional) a key-value pair of an edge attribute and a predicate function. Bonds with attributes that cause the predicate to return false are excluded from writing.
"""
function write_bond_information(crystal::Crystal, filename::String;
        center_at_origin::Bool=false, bond_filter::Pair{Symbol, F}=(:NOTHING => x -> ())) where F <: Function
    if ne(crystal.bonds) == 0
        @warn("Crystal %s has no bonds present. To get bonding information for this
        crystal run `infer_bonds!` with an array of bonding rules\n", crystal.name)
    end
    if ! occursin(".vtk", filename)
        filename *= ".vtk"
    end
    # filter bonds
    idx_keep_bonds = trues(ne(crystal.bonds))
    attr, pred = bond_filter
    if attr != :NOTHING
        for (b, bond) ∈ enumerate(edges(crystal.bonds))
            prop = get_prop(crystal.bonds, bond, attr)
            if !pred(prop)
                idx_keep_bonds[b] = false
            end
        end
    end
    # write output
    vtk_file = open(filename, "w")
    @printf(vtk_file, "# vtk DataFile Version 2.0\n%s bond information\nASCII\n
        DATASET POLYDATA\nPOINTS %d double\n", crystal.name, nv(crystal.bonds))
    for i = 1:crystal.atoms.n
        if center_at_origin
            @printf(vtk_file, "%0.5f\t%0.5f\t%0.5f\n", (crystal.box.f_to_c *
                (crystal.atoms.coords.xf[:, i] - [0.5, 0.5, 0.5]))...)
        else
            @printf(vtk_file, "%0.5f\t%0.5f\t%0.5f\n", (crystal.box.f_to_c *
                crystal.atoms.coords.xf[:, i])...)
        end
    end
    @printf(vtk_file, "\nLINES %d %d\n", sum(idx_keep_bonds), 3 * sum(idx_keep_bonds))
    for (e, edge) in enumerate(edges(crystal.bonds))
        if idx_keep_bonds[e]
            @printf(vtk_file, "2\t%d\t%d\n", edge.src - 1, edge.dst - 1)
        end
    end
    close(vtk_file)
    @printf("Saving bond information for crystal %s to %s.\n", crystal.name,
        joinpath(pwd(), filename))
end


"""
Loop through xtal and calculate any missing distances
"""
function calc_missing_bond_distances!(xtal::Crystal)
    for bond in collect(edges(xtal.bonds))
        if !(:distance ∈ keys(props(xtal.bonds, bond)))
            i = src(bond)
            j = dst(bond)
            set_prop!(xtal.bonds, i, j, :distance, distance(xtal.atoms, xtal.box, i, j, true))
        end
    end
end


"""
Gets rid of the bonds across unit cell boundaries
"""
function drop_cross_pb_bonds!(bonds::MetaGraph)
    for bond in collect(edges(bonds))
        if get_prop(bonds, bond, :cross_boundary)
            rem_edge!(bonds, src(bond), dst(bond))
        end
    end
end

drop_cross_pb_bonds!(xtal::Crystal) = drop_cross_pb_bonds!(xtal.bonds)


"""
    `clear_vectors!(xtal)`
    `clear_vectors!(xtal.bonds)`

Removes edge vector attributes from crystal bonding graph.
"""
function clear_vectors!(bonds::MetaGraph)
    for edge ∈ edges(bonds)
        rem_prop!(bonds, edge, :vector)
        rem_prop!(bonds, edge, :direction)
    end
    set_prop!(bonds, :has_vectors, false)
end

clear_vectors!(xtal::Crystal) = clear_vectors!(xtal.bonds)


"""
    `calculate_bond_vectors!(xtal)`
Adds a property to the edges of a crystal's bonding graph giving the vector between source and destination nodes in Cartesian space.
"""
function calculate_bond_vectors!(xtal::Crystal)
    # check for bonds
    if ne(xtal.bonds) == 0
        # no bonds--but could be legit property of xtal. non-fatal error
        @error "Crystal $(xtal.name) has no bonds."
        return
    end
    # check that vectors have not already been calculated
    if has_prop(xtal.bonds, :has_vectors) && get_prop(xtal.bonds, :has_vectors)
        # loop will skip all bonds if vectors already set. fatal error.
        error("Crystal $(xtal.name) already has vectors. Clear them with `clear_vectors!` first.")
        return
    end
    # get bonds as directed graph (matters for vectors!)
    bonds = DiGraph(adjacency_matrix(xtal.bonds))
    # loop over bonding edges
    for edge ∈ edges(bonds)
        # get node IDs
        i, j = src(edge), dst(edge)
        # only need to label an edge once in undirected graph
        if has_prop(xtal.bonds, j, i, :vector)
            continue
        end
        # get fractional coordinates, vector
        xf_i = xtal.atoms.coords.xf[:, i]
        xf_j = xtal.atoms.coords.xf[:, j]
        dxf = xf_j .- xf_i
        # apply nearest image convention to dxf
        nearest_image!(dxf)
        # transform to Cartesian
        dxc = xtal.box.f_to_c * dxf
        # annotate edge in graph
        set_prop!(xtal.bonds, i, j, :vector, dxc)
        # xtal.bonds is undirected, so need to specify direction of vector
        set_prop!(xtal.bonds, i, j, :direction, (i, j))
    end
    set_prop!(xtal.bonds, :has_vectors, true)
end


"""
    `get_bond_vector(bonds, i, j)`
Returns the vector between connected nodes `i` and `j` of the `bonds` graph.
"""
function get_bond_vector(bonds::MetaGraph, i::Int, j::Int)
    if get_prop(bonds, i, j, :direction) == (i, j)
        return get_prop(bonds, i, j, :vector)
    else
        return -1 .* get_prop(bonds, i, j, :vector)
    end
end

get_bond_vector(xtal::Crystal, i::Int, j::Int) = get_bond_vector(xtal.bonds, i, j)


"""
    `bond_angle(xtal.bonds, 2, 1, 3)`
    `bond_angle(xtal, 8, 121, 42)`

Returns the bond angle between three nodes of a bonding graph (or three atoms in a crystal), if the edges (bonds) exist.
Otherwise, returns `NaN`
"""
function bond_angle(bonds::MetaGraph, i::Int, j::Int, k::Int)::Float64
    # θ = arccos(u⃗•v⃗ / (‖u⃗‖*‖v⃗‖))
    @assert has_prop(bonds, :has_vectors) "Vectors have not been calculated. Do so with `calculate_bond_vectors!`"
    @assert has_edge(bonds, i, j) && has_edge(bonds, j, k) "Invalid vertex selection. Cannot calculate ∠($i,$j,$k)"
    # get vectors in correct orientation
    u⃗ = get_bond_vector(bonds, j, i)
    v⃗ = get_bond_vector(bonds, j, k)
    norm_u⃗ = get_prop(bonds, i, j, :distance)
    norm_v⃗ = get_prop(bonds, j, k, :distance)
    @assert norm_u⃗ ≠ 0 && norm_v⃗ ≠ 0
    @assert !isnan(norm_u⃗) && !isnan(norm_v⃗)
    ϕ = dot(u⃗, v⃗) / (norm_u⃗ * norm_v⃗)
    # correct for numerical precision domain error |ϕ| > 1
    if isapprox(ϕ, 1)
        ϕ = 1.0
    elseif isapprox(ϕ, -1)
        ϕ = -1.0
    end
    return acos(ϕ)
end

bond_angle(xtal::Crystal, i::Int, j::Int, k::Int) = bond_angle(xtal.bonds, i, j, k)


"""
    `bond_distance(xtal, i, j)`
    `bond_distance(xtal.bonds, i, j)`

Gives the distance between two bonded atoms in a crystal or two nodes in a bonding graph.
Returns `NaN` if the atoms (nodes) are not bonded (connected).
"""
function bond_distance(bonds::MetaGraph, i::Int, j::Int)::Float64
    if has_edge(bonds, i, j)
        return get_prop(bonds, i, j, :distance)
    else
        return NaN
    end
end

bond_distance(xtal::Crystal, i::Int, j::Int) = bond_distance(xtal.bonds, i, j)
