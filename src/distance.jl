"""
    nearest_image!(dxf)

Applies the nearest image convention on a vector `dxf` between two atoms in fractional
space; modifies `dxf` for nearest image convention. Fractional coordinates here fall in
[0, 1] so that the box is [0, 1]^3 in fractional space.

Warning: this assumes the two molecules are in the box described by fractional coords [0, 1]³.

# Arguments
- `dxf::Array{Float64}`: A vector between two atoms in fractional space
"""
@inline function nearest_image!(dxf::Array{Float64})
    for i in eachindex(dxf)
        @inbounds if abs(dxf[i]) > 0.5
            @inbounds dxf[i] -= sign(dxf[i])
        end
    end
end

@doc raw"""
    r = distance(coords, box, i, j, apply_pbc)
    r = distance(atoms, box, i, j, apply_pbc) # atoms i and j
    r = distance(charges, box, i, j, apply_pbc) # charges i and j
    r = distance(atoms, i, j) # no PBCs, coords must be in Cartesian coords
    r = distance(coords, i, j) # no PBCs, coords must be in Cartesian coords

calculate the (Cartesian) distance between particles `i` and `j`.

apply periodic boundary conditions if and only if `apply_pbc` is `true`.

# arguments
- `coords::Coords`: the coordinates (`Frac>:Coords` or `Cart>:Coords`)
- `atoms::Atoms`: atoms
- `charges::charges`: atoms
- `box::Box`: unit cell information
- `i::Int`: index of the first particle
- `j::Int`: Index of the second particle
- `apply_pbc::Bool`: `true` if we wish to apply periodic boundary conditions, `false` otherwise
"""
function distance(coords::Frac, box::Box, i::Int, j::Int, apply_pbc::Bool)
    dxf = coords.xf[:, i] - coords.xf[:, j]
    if apply_pbc
        nearest_image!(dxf)
    end
    return norm(box.f_to_c * dxf)
end

function distance(coords::Cart, box::Box, i::Int, j::Int, apply_pbc::Bool)
    dx = coords.x[:, i] - coords.x[:, j]
    if apply_pbc
        dxf = box.c_to_f * dx
        nearest_image!(dxf)
        return norm(box.f_to_c * dxf)
    else
        return norm(dx)
    end
end

# no PBCs
function distance(coords::Cart, i::Int, j::Int)
    dx = coords.x[:, i] - coords.x[:, j]
    return norm(dx)
end
distance(atoms::Atoms{Cart}, i::Int, j::Int) = distance(atoms.coords, i, j)

distance(atoms::Atoms, box::Box, i::Int, j::Int, apply_pbc::Bool) = distance(atoms.coords, box, i, j, apply_pbc)
distance(charges::Charges, box::Box, i::Int, j::Int, apply_pbc::Bool) = distance(charges.coords, box, i, j, apply_pbc)

function pairwise_distances(coords::Frac, box::Box, apply_pbc::Bool)
    n = length(coords)
    pd = zeros(n, n)
    for i = 1:n
        for j = 1:n
            pd[i, j] = distance(coords, box, i, j, apply_pbc)
            if i > j
                pd[i, j] =  pd[j, i]
            end
        end
    end
    return pd
end
pairwise_distances(coords::Cart, box::Box, apply_pbc::Bool) = pairwise_distances(Frac(coords, box), box, apply_pbc)

"""
    overlap_flag, overlap_pairs = overlap(crystal; r_crit=nothing, apply_pbc=true)

determine the pairs of atoms (if any) in a crystal that overlap.

two atoms are defined to overlap if their (Cartesian) distance is less than `r_crit` Å.

if `isnothing(r_crit)`, we use covalent radii to determine overlap.
covalent radii are obtained from [`get_covalent_radii`](@ref).
if covalent radii not available for an atom, we use 0.3 Å

# Arguments
- `xtal::Crystal`: the crystal
- `r_crit::Union{Float64, Nothing}`: atoms overlap if a distances less than `r_crit` apart. if `isnothing(r_crit)`, use covalent radii.
- `apply_pbc::Bool`: `true` if we wish to apply periodic boundary conditions, `false` otherwise

# Returns
* `overlap_flag::Bool`: `true` if overlap, `false` otherwise
* `overlap_ids::Array{Tuple{Int, Int}, 1}`: ids of atom pairs that are overlapping e.g. `[(4, 5), (7, 8)]`
"""
function overlap(xtal::Crystal; r_crit::Union{Float64, Nothing}=nothing, apply_pbc::Bool=true, throw_error_print_warnings::Bool=true)
    overlap_flag = false
    overlap_ids = Array{Tuple{Int, Int}, 1}()
    covalent_radii = get_covalent_radii()
    
    # loop over pairs of atoms in the crystal.
    for i = 1:xtal.atoms.n
        for j = i+1:xtal.atoms.n
            atoms_ij_overlap = false

            # compute distance between the pair of atoms
            r = distance(xtal.atoms.coords, xtal.box, i, j, apply_pbc)

            if isnothing(r_crit)
                # use covalent radii if available
                species_i = strip_number_from_label(xtal.atoms.species[i])
                species_j = strip_number_from_label(xtal.atoms.species[j])
                
                # use smallest radii to avoid false overlaps
                r_i = 0.3
                r_j = 0.3
                if species_i in keys(covalent_radii)
                    r_i = covalent_radii[species_i][:radius_Å] - covalent_radii[species_i][:esd_Å]
                end
                if species_j in keys(covalent_radii)
                    r_j = covalent_radii[species_j][:radius_Å] - covalent_radii[species_j][:esd_Å]
                end
                if r < r_i + r_j - 0.25
                    atoms_ij_overlap = true
                end
            else
                # use provided r_crit
                if r < r_crit
                    atoms_ij_overlap = true
                end
            end
            
            if atoms_ij_overlap
                push!(overlap_ids, (i, j))
                overlap_flag = true
                if throw_error_print_warnings
                    @warn @sprintf("atom %d (%s) and %d (%s) in %s are overlapping. r = %f Å\n", i, xtal.atoms.species[i], 
                                                                                       j, xtal.atoms.species[j],
                                                                                       xtal.name, r)
                end
            end
        end
    end
    if overlap_flag && throw_error_print_warnings
        error(xtal.name * " has overlapping pairs of atoms!
            (this occurs often when applying symmetry rules, if your structure was not in P1 symmetry to begin with.)
            to investigate or avoid this error, pass `check_overlap=false` to the `Crystal` constructor, 
            then run `overlap(crystal, throw_error_print_warnings=false)` to find pairs of overlapping atoms for inspection.
            or pass `remove_duplicates=true` to the `Crystal` constructor to automatically remove the duplicate atoms and charges.
            you should then visualize your .cif to make sure it is proper.`\n")
    end
    return overlap_flag, overlap_ids
end

"""
    atoms = remove_duplicates(atoms, box, apply_pbc, r_tol=0.1)
    charges = remove_duplicates(charges, box, apply_pbc, r_tol=0.1)

remove duplicates from atoms or charges.

loops through all pairs of atoms/charges. if a pair is a duplicate, one is deleted.

two atoms are duplicates if both:
* same species
* less than a distance `r_tol` apart
two charges are duplicate if both:
* charge values are within `q_tol`
* less than a distance `r_tol` apart

# arguments
- `atoms::Atoms`: the atoms
- `charges::Charges`: the charges
- `box::Box`: unit cell information
- `apply_pbc::Bool`: true iff we apply periodic boundary conditions when computing the distance.
- `r_tol::Float64`: atoms/charges are overlapping if within `r_tol` distance (PBC applied)
- `q_tol::Float64`: charges have the same charge value if their charges are within `q_tol` of each other
"""
function remove_duplicates(ac::Union{Atoms{Frac}, Charges{Frac}}, box::Box, apply_pbc::Bool;
                           r_tol::Float64=0.1, q_tol::Float64=0.0001)
    ids_keep = trues(ac.n)
    for i = 1:ac.n
        for j = i+1:ac.n
            if isa(ac, Atoms{Frac})
                # if different species, not duplicates.
                if ac.species[i] != ac.species[j]
                    continue
                end
            elseif isa(ac, Charges{Frac})
                # if different value of charge, not duplicates.
                if ! isapprox(ac.q[i], ac.q[j], atol=q_tol)
                    continue
                end
            end
            # are they overlapping?
            r = distance(ac.coords, box, i, j, apply_pbc)
            if r < r_tol
                ids_keep[j] = false
            end
        end
    end
    return ac[ids_keep]
end


"""
    dm = distance_matrix(crystal, apply_pbc)

Compute the distance matrix `a` of the crystal, where `a[i, j]` is the
distance between atom `i` and `j`. This matrix is symmetric. If `apply_pbc = true`,
periodic boundary conditions are applied when computing the distance.

# Arguments
-`crystal::Crystal`: crystal structure
-`apply_pbc::Bool`: whether or not to apply periodic boundary conditions when computing the distance

# Returns
-`dm::Array{Float64, 2}`: symmetric, square distance matrix with zeros on the diagonal
"""
function distance_matrix(crystal::Crystal, apply_pbc::Bool)
    dm = zeros(crystal.atoms.n, crystal.atoms.n)
    for i = 1:crystal.atoms.n
        for j = (i+1):crystal.atoms.n
            dm[i, j] = distance(crystal.atoms, crystal.box, i, j, apply_pbc)
            dm[j, i] = dm[i, j] # symmetric
        end
    end
    return dm
end
