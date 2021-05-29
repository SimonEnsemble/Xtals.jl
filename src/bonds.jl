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
function bondingrules(;covalent_radii::Dict{Symbol,Float64}=get_global(:covalent_radii), pad::Float64=0.)::Array{BondingRule}
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
            push!(bondingrules, BondingRule(atom_i, atom_j, max_dist))
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
function write_bonding_rules(filename::String, bonding_rules::Union{Array{BondingRule},Nothing}=nothing)
    bonding_rules = isnothing(bonding_rules) ? get_global(:bonding_rules) : bonding_rules
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
    br = get_global(:bonding_rules)
    set_global(:bonding_rules => vcat(bonding_rules, br))
end
add_bonding_rules(bonding_rule::BondingRule) = add_bonding_rules([bonding_rule])


# for pretty-printing the bonding rules
function show(io::IO, bonding_rules::Array{BondingRule})
    for r in bonding_rules
        println("%s\t%s\t%.3f\t%.3f", r.species_i, r.species_j, r.max_dist)
    end
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
                return true, r, cross_boundary_flag
            else
                # return (don't continue!) b/c we found relevant bonding rule, don't apply others
                return false, NaN, cross_boundary_flag
            end
        end
    end
    # no bonding rule was applied; therefore not bonded.
    @warn "No bonding rule for $(species_i)-$(species_j) interaction."
    return false, NaN, cross_boundary_flag
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
end


"""
    infer_bonds!(crystal, include_bonds_across_periodic_boundaries, bonding_rules=nothing)

Populate the bonds in the crystal object based on the bonding rules. If a
pair doesn't have a suitable rule then they will not be considered bonded.

`:*` is considered a wildcard and can be substituted for any species. It is a
good idea to include a bonding rule between two `:*` to allow any atoms to bond
as long as they are close enough.

The bonding rules are hierarchical, i.e. the first bonding rule takes precedence over the latter ones.

# Arguments
- `crystal::Crystal`: The crystal that bonds will be added to
- `include_bonds_across_periodic_boundaries::Bool`: Whether to check across the periodic boundary when calculating bonds
- `bonding_rules::Union{Array{BondingRule, 1}, Nothing}=nothing`: The array of bonding rules that will be used to fill the bonding information. They are applied in the order that they appear. if `nothing`, default bonding rules will be applied; see [`get_bonding_rules`](@ref)
"""
function infer_bonds!(crystal::Crystal, include_bonds_across_periodic_boundaries::Bool;
                      bonding_rules::Union{Array{BondingRule, 1}, Nothing}=nothing)
    @assert ne(crystal.bonds) == 0 @sprintf("The crystal %s already has bonds. Remove them with the `remove_bonds!` function before inferring new ones.", crystal.name)
    bonding_rules = isnothing(bonding_rules) ? get_global(:bonding_rules) : bonding_rules
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
    bond_sanity_check(crystal)
end


"""
    ids_neighbors, xs, rs = neighborhood(crystal, i, r, dm)

Find and characterize the neighborhood of atom `i` in the crystal `crystal`.
A neighborhood is defined as all atoms within a distance `r` from atom `i`.
The distance matrix `dm` is used to find the distances of all other atoms in the crystal from atom `i`.

# Arguments
- `crystal::Crystal`: crystal structure
- `i::Int`: Index of the atom (in `crystal`) which the neighborhood is to be characterized.
- `r::Float64`: The maximum distance the neighborhood will be characterized.
- `dm::Array{Float64, 2}`: The distance matrix, see [`distance_matrix`](@ref)

# Returns
- `ids_neighbors::Array{Int, 1}`: indices of `crystal.atoms` within the neighborhood of atom `i`.
- `xs::Array{Array{Float64, 1}, 1}`: array of Cartesian positions of the atoms surrounding atom `i`.
    The nearest image convention has been applied to find the nearest periodic image. Also, the coordinates of atom `i`
    have been subtracted off from these coordinates so that atom `i` lies at the origin of this new coordinate system.
    The first vector in `xs` is `[0, 0, 0]` corresponding to atom `i`.
    The choice of type is for the Voronoi decomposition in Scipy.
- `rs::Array{Float64, 1}`: list of distances of the neighboring atoms from atom `i`.
"""
function neighborhood(crystal::Crystal, i::Int, r::Float64, dm::Array{Float64, 2})
    # get indices of atoms within a distance r of atom i
    #  the greater than zero part is to not include itself
    ids_neighbors = findall((dm[:, i] .> 0.0) .& (dm[:, i] .< r))
    # rs is the list of distance of these neighbors from atom i
    rs = [dm[i, id_n] for id_n in ids_neighbors]
    @assert all(rs .< r)
    # xs is a list of Cartesian coords of the neighborhood
    #  coords of atom i are subtracted off
    #  first entry is coords of atom i, the center, the zero vector
    #  remaining entries are neighbors
    # this list is useful to pass to Voronoi for getting Voronoi faces
    #  of the neighborhood.
    xs = [[0.0, 0.0, 0.0]] # this way atom zero is itself
    for j in ids_neighbors
        # subtract off atom i, apply nearest image
        xf = crystal.atoms.coords.xf[:, j] - crystal.atoms.coords.xf[:, i]
        nearest_image!(xf)
        x = crystal.box.f_to_c * xf
        push!(xs, x)
    end
    return ids_neighbors, xs , rs
end


"""
    ids_shared_voro_face = _shared_voronoi_faces(ids_neighbors, xs)

Of the neighboring atoms, find those that share a Voronoi face.

# Arguments
- `ids_neighbors::Array{Int, 1}`: indices of atoms within the neighborhood of a specific atom.
- `xs::Array{Array{Float64, 1}, 1}`: array of Cartesian position of the atoms within the neighborhood of a specific atom, relative to the specific atom.

# Returns
- `ids_shared_voro_face::Array{Int, 1}`: indices of atoms that share a Voronoi face with a specific atom
"""
function _shared_voronoi_faces(ids_neighbors::Array{Int,1}, xs::Array{Array{Float64,1},1})
    scipy = pyimport("scipy.spatial")
    # first element of xs is the point itself, the origin
    @assert length(ids_neighbors) == (length(xs) - 1)
    voro = scipy.Voronoi(xs)
    rps = voro.ridge_points # connection with atom zero are connection with atom i
    ids_shared_voro_face = Int[] # corresponds to xs, not to atoms of crystal
    for k = 1:size(rps)[1]
        if sort(rps[k, :])[1] == 0 # a shared face with atom i!
            push!(ids_shared_voro_face, sort(rps[k, :])[2])
        end
    end
    # zero based indexing in Scipy accounted for since xs[0] is origin, atom i.
    return ids_neighbors[ids_shared_voro_face]
end


"""
    ids_bonded = bonded_atoms(crystal, i, dm; r=6., σ=3.)

Returns the ids of atoms that are bonded to atom `i` by determining bonds using a Voronoi method and covalent radius data (see [`covalent_radii`](@ref))

# Arguments
- `crystal::Crystal`: Crystal structure in which the bonded atoms will be determined
- `i::Int`: Index of the atom we want to determine the bonds of
- `dm::Array{Float64, 2}`: The distance matrix, see [`distance_matrix`](@ref)
- `r::Float64`: The maximum distance used to determine the neighborhood of atom `i`
- `covalent_radii::Dict{Symbol, Float64}`: Cordero parameter dictionary. See [`covalent_radii`](@ref)
- `pad::Float64`: The amount to pad the covalent radius in Å

# Returns
- `ids_bonded::Array{Int, 1}`: A list of indices of atoms bonded to atom `i`
"""
function bonded_atoms(crystal::Crystal, i::Int, dm::Array{Float64, 2},
        r::Float64, pad::Float64, covalent_radii::Dict{Symbol, Float64})
    species_i = crystal.atoms.species[i]
    ids_neighbors, xs, _ = neighborhood(crystal, i, r, dm)
    ids_shared_voro_faces = _shared_voronoi_faces(ids_neighbors, xs)
    ids_bonded = Int[]
    for j in ids_shared_voro_faces
        species_j = crystal.atoms.species[j]
        # sum of covalent radii
        radii_sum = covalent_radii[species_j] + covalent_radii[species_i]
        # margin = σ e.s.d.s, unless that's too small
        max_dist = radii_sum + pad
        if dm[i, j] ≤ max_dist
            push!(ids_bonded, j)
        end
    end
    return ids_bonded
end


"""
    infer_geometry_based_bonds!(crystal, include_bonds_across_periodic_boundaries::Bool)

Infers bonds by first finding which atoms share a Voronoi face, and then bond the atoms if the distance
 between them is less than the sum of the covalent radius of the two atoms (plus a tolerance).

# Arguments
- `crystal::Crystal`: The crystal structure
- `include_bonds_across_periodic_boundaries::Bool`: Whether to check across the periodic boundaries
- `r::Float`: voronoi radius, Å
- `pad::Float`: amount to pad covalent radii, Å
- `covalent_radii::Dict{Symbol,Float64}`: See [`covalent_radii`](@ref)
"""
function infer_geometry_based_bonds!(crystal::Crystal,
        include_bonds_across_periodic_boundaries::Bool;
        r::Float64=6., pad::Float64=0.,
        covalent_radii::Union{Nothing,Dict{Symbol,Float64}}=nothing)
    @assert ne(crystal.bonds) == 0 @sprintf("The crystal %s already has bonds. Remove them with the `remove_bonds!` function before inferring new ones.", crystal.name)
    if isnothing(covalent_radii)
        covalent_radii = get_global(:covalent_radii)
    end
    dm = distance_matrix(crystal, include_bonds_across_periodic_boundaries)
    for i = 1:crystal.atoms.n
        for j in bonded_atoms(crystal, i, dm, r, pad, covalent_radii)
            add_edge!(crystal.bonds, i, j)
        end
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
    write_bond_information(xtal, filename, :cross_boundary => p -> p, "bonds.vtk") # cross boundary bonds only
    write_bond_information(xtal, filename, :distance => d -> d < 1.0, "bonds.vtk") # distance less than 1.0

Writes the bond information from a crystal to the selected filename.

# Arguments
- `crystal::Crystal`: The crystal to have its bonds written to a vtk file
- `filename::String`: The filename the bond information will be saved to. If left out, will default to crystal name.
- `center_at_origin::Bool`: (optional) center the coordinates at the origin of the crystal
- `bond_filter::Pair{Symbol, Function}`: (optional) a key-value pair of an edge attribute and a predicate function. Bonds with attributes that cause the predicate to return false are excluded from writing.
"""
function write_bond_information(crystal::Crystal, filename::String;
        center_at_origin::Bool=false, bond_filter::Union{Pair{Symbol, F}, Nothing}=nothing) where F
    if ne(crystal.bonds) == 0
        @warn("Crystal %s has no bonds present. To get bonding information for this
        crystal run `infer_bonds!` with an array of bonding rules\n", crystal.name)
    end
    if ! occursin(".vtk", filename)
        filename *= ".vtk"
    end
    # filter bonds
    idx_keep_bonds = trues(ne(crystal.bonds))
    if !isnothing(bond_filter)
        attr, pred = bond_filter
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
        p = get_prop(bonds, bond, :cross_boundary)
        if ismissing(p)
            @warn "Bond w/ missing :cross_boundary attribute"
        elseif p
            rem_edge!(bonds, src(bond), dst(bond))
        end
    end
end

drop_cross_pb_bonds!(xtal::Crystal) = drop_cross_pb_bonds!(xtal.bonds)
