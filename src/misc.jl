function add_extension(filename::String, extension::String)
    if ! occursin(extension, filename)
        filename *= extension
    end
    return filename
end

"""
    atoms = read_xyz("molecule.xyz")

read a list of atomic species and their corresponding coordinates from an .xyz file.

# Arguments
- `filename::AbstractString`: the path to and filename of the .xyz file

# Returns
- `atoms::Atoms{Cart}`: the set of atoms read from the .xyz file.
"""
function read_xyz(filename::AbstractString)
    f = open(filename)
    lines = readlines(f)
    if length(lines) == 0
        return Symbol[], Float64[]
    end
    n = parse(Int, lines[1]) # get number of atoms
    species = Symbol[]
    x = zeros(Float64, 3, n)
    for i = 1:n
        push!(species, Symbol(split(lines[i + 2])[1]))
        for j = 1:3
            x[j, i] = parse(Float64, split(lines[i + 2])[1 + j])
        end
    end
    close(f)
    return Atoms(species, Cart(x))
end

"""
    write_xyz(atoms, filename; comment="")
    write_xyz(crystal; comment="", center_at_origin=false)
    write_xyz(molecules, box, filename; comment="") # fractional
    write_xyz(molecules, box, filename; comment="") # Cartesian

write atoms to an .xyz file.

# Arguments
- `atoms::Atoms`: the set of atoms.
- `filename::AbstractString`: the filename (absolute path) of the .xyz file. (".xyz" appended automatically
if the extension is not provided.)
- `comment::AbstractString`: comment if you'd like to write to the file.
- `center_at_origin::Bool`: (for crystal only) if `true`, translate all coords such that the origin is the center of the unit cell.
"""
function write_xyz(atoms::Atoms{Cart}, filename::AbstractString; comment::AbstractString="")
    filename = add_extension(filename, ".xyz")

    xyzfile = open(filename, "w")
    @printf(xyzfile, "%d\n%s\n", atoms.n, comment)
    for i = 1:atoms.n
		@printf(xyzfile, "%s\t%.4f\t%.4f\t%.4f\n", atoms.species[i],
            atoms.coords.x[1, i], atoms.coords.x[2, i], atoms.coords.x[3, i])
    end
    close(xyzfile)
    return nothing
end


"""
    atoms, bonds = read_mol("molecule.mol")

read a `.mol` file, which contains info about both atoms and bonds.
see [here](https://chem.libretexts.org/Courses/University_of_Arkansas_Little_Rock/ChemInformatics_(2017)%3A_Chem_4399%2F%2F5399/2.2%3A_Chemical_Representations_on_Computer%3A_Part_II/2.2.2%3A_Anatomy_of_a_MOL_file) for the anatomy of a `.mol` file.

# Arguments
- `filename::AbstractString`: the path to and filename of the .mol file (must pass extension)

# Returns
- `atoms::Atoms{Cart}`: the set of atoms read from the `.mol` file.
- `bonds::MetaGraph`: the bonding graph of the atoms read from the `.mol` file.
- `bond_types::Array{Int, 1}`: the array of bond types.
"""
function read_mol(filename::String)
    f = open(filename)
    lines = readlines(f)
    close(f)

    n_atoms = parse(Int, split(lines[4])[1])
    n_bonds = parse(Int, split(lines[4])[2])

    atoms = Atoms([:blah for i = 1:n_atoms],
                  Cart(zeros(Float64, 3, n_atoms))
                  )
    for i = 1:n_atoms
        line_atom_i = split(lines[4+i])

        atoms.species[i] = Symbol(line_atom_i[4])
        for j = 1:3
            atoms.coords.x[j, i] = parse(Float64, line_atom_i[j])
        end
    end

    bonds = MetaGraph(n_atoms)
    bond_types = [-1 for i = 1:n_bonds]
    for b = 1:n_bonds
        line_bond_b = split(lines[n_atoms + 4 + b])

        i = parse(Int, line_bond_b[1])
        j = parse(Int, line_bond_b[2])
        add_edge!(bonds, i, j)

        bond_types[b] = parse(Int, line_bond_b[3])
    end

    return atoms, bonds, bond_types
end
