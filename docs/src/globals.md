# Global Variables

In `Xtals.jl`, global variables are stored in a dictionary called `rc`.  The entries in `rc` are used for holding various important pieces of information:

- `rc[:atomic_masses]` : a dictionary which maps atomic species to their masses in AMU
- `rc[:cpk_colors]` : a dictionary which maps atomic species to their CPK colors
- `rc[:bonding_rules]` : an array which lists the maximum bonding distances between each pair of atom types in Å
- `rc[:covalent_radii]` : a dictionary which maps atomic species to their covalent radii in Å
- `rc[:paths]` : a dictionary of paths to directories to search for input data

The `rc[:paths]` dictionary holds the following:

- `rc[:paths][:data]` : not directly used by `Xtals.jl`, but important for other packages which use `Xtals`, e.g. [PorousMaterials.jl](https://simonensemble.github.io/PorousMaterials.jl/dev/)
- `rc[:paths][:crystals]` : where to find crystallographic input data
