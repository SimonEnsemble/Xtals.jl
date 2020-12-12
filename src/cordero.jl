"""
    covalent_radii = get_covalent_radii()

Create a dictionary with the Cordero covalent radius and estimated standard deviation for each element, using the data in `PATH_TO_DATA/cordero.csv`

    covalent_radii = get_covalent_radii("my_params.csv")

Create a dictionary with the Cordero covalent radius and estimated standard deviation for each element specified in `PATH_TO_DATA/my_params.csv`

# Arguments
-`cordero_data::String`: name of file containing covalent radii and estimated standard deviations.

# Returns
-`covalent_radii::Dict{Symbol, Dict{Symbol, Float64}}`: A dictionary with elements as keys and dictionaries of respective cordero covalent radii and e.s.d.s as the values.

# Example
```julia
covalent_radii = get_covalent_radii()
covalent_radii[:N][:radius_Å] # 0.71
covalent_radii[:N][:esd_pm] # 1.0
```
"""
function get_covalent_radii(; cordero_data::Union{String,Nothing}=nothing)
    if cordero_data == nothing
        df = COVALENT_RADII
    else
        # read Cordero data
        df = CSV.read(joinpath(PATH_TO_DATA, cordero_data), DataFrame, comment="#")
    end
    # parse into params dict
    covalent_radii = Dict{Symbol, Dict{Symbol, Float64}}()
    for atom in eachrow(df)
        covalent_radii[Symbol(atom[:atom])] = Dict(:radius_Å => atom.covalent_radius_A, :esd_pm => atom.esd_pm)
    end
    # Carbon, Iron, Manganese, and Cobalt have multiple entries due to hybridization/spin
    # Generic Carbon approximation; sp3 hybridization is the largest r, sp2 is the largest esd. Mix-n-match.
    covalent_radii[:C] = Dict(:radius_Å=> covalent_radii[:C_sp3][:radius_Å], :esd_pm => covalent_radii[:C_sp2][:esd_pm])
    # Generic Mn, Fe, Co approx's. High-spin r and esd is larger for each. Use high-spin.
    for x in [:Mn, :Fe, :Co]
        covalent_radii[x] = covalent_radii[Symbol(string(x) * "_hi")]
    end
    return covalent_radii
end


# makes dataframe out of stringified cordero.csv
function _CorderoParams(headers::Array{Symbol}, data::String)::DataFrame
    species = Symbol[]
    radius = Float64[]
    esd = Float64[]
    for line ∈ split(data,"\n")
        fields = split(line, ",")
        if length(fields) == 3
            push!(species, Symbol(fields[1]))
            push!(radius, parse(Float64, fields[2]))
            push!(esd, parse(Float64, fields[3]))
        end
    end
    return DataFrame([species, radius, esd], headers)
end

# Cordero bond radii, DOI: 10.1039/B801115J Table 2.
COVALENT_RADII = _CorderoParams([:atom, :covalent_radius_A, :esd_pm],
"""H,0.31,5
He,0.28,1
Li,1.28,7
Be,0.96,3
B,0.84,3
C_sp3,0.76,1
C_sp2,0.73,2
C_sp,0.69,1
N,0.71,1
O,0.66,2
F,0.57,3
Ne,0.58,1
Na,1.66,9
Mg,1.41,7
Al,1.21,4
Si,1.11,2
P,1.07,3
S,1.05,3
Cl,1.02,4
Ar,1.06,10
K,2.03,12
Ca,1.76,10
Sc,1.7,7
Ti,1.6,8
V,1.53,8
Cr,1.39,5
Mn_lo,1.39,5
Mn_hi,1.61,8
Fe_lo,1.32,3
Fe_hi,1.52,6
Co_lo,1.26,3
Co_hi,1.5,7
Ni,1.24,4
Cu,1.32,4
Zn,1.22,4
Ga,1.22,3
Ge,1.2,4
As,1.19,4
Se,1.2,4
Br,1.2,3
Kr,1.16,4
Rb,2.2,9
Sr,1.95,10
Y,1.9,7
Zr,1.75,7
Nb,1.64,6
Mo,1.54,5
Tc,1.47,7
Ru,1.46,7
Rh,1.42,7
Pd,1.39,6
Ag,1.45,5
Cd,1.44,9
In,1.42,5
Sn,1.39,4
Sb,1.39,5
Te,1.38,4
I,1.39,3
Xe,1.4,9
Cs,2.44,11
Ba,2.15,11
La,2.07,8
Ce,2.04,9
Pr,2.03,7
Nd,2.01,6
Pm,1.99,1
Sm,1.98,8
Eu,1.98,6
Gd,1.96,6
Tb,1.94,5
Dy,1.92,7
Ho,1.92,7
Er,1.89,6
Tm,1.9,10
Yb,1.87,8
Lu,1.87,8
Hf,1.75,10
Ta,1.7,8
W,1.62,7
Re,1.51,7
Os,1.44,4
Ir,1.41,6
Pt,1.36,5
Au,1.36,6
Hg,1.32,5
Tl,1.45,7
Pb,1.46,5
Bi,1.48,4
Po,1.4,4
At,1.5,1
Rn,1.5,1
Fr,2.6,1
Ra,2.21,2
Ac,2.15,1
Th,2.06,6
Pa,2,1
U,1.96,7
Np,1.9,1
Pu,1.87,1
Am,1.8,6
Cm,1.69,3
""")
