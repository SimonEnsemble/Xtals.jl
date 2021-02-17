"""
    covalent_radii = get_covalent_radii()

Create a dictionary with the Cordero covalent radius and estimated standard deviation for each element, using the data in `PATH_TO_DATA/cordero.csv`

    covalent_radii = get_covalent_radii()

Create a dictionary with the Cordero covalent radius and estimated standard deviation for each element specified in `PATH_TO_DATA/my_params.csv`

# Returns
- `covalent_radii::Dict{Symbol, Dict{Symbol, Float64}}`: A dictionary with elements as keys and dictionaries of respective cordero covalent radii and e.s.d.s as the values.

# Example
```julia
covalent_radii = get_covalent_radii()
covalent_radii[:N][:radius_Å] # 0.71
covalent_radii[:N][:esd_Å] # 0.01
```
"""
function get_covalent_radii()
    # parse into params dict
    covalent_radii = Dict{Symbol, Dict{Symbol, Float64}}()
    for atom in eachrow(COVALENT_RADII)
        covalent_radii[Symbol(atom[:atom])] = Dict(:radius_Å => atom.covalent_radius_A, :esd_Å => atom.esd_Å)
    end
    # Carbon, Iron, Manganese, and Cobalt have multiple entries due to hybridization/spin
    # Generic Carbon approximation; sp3 hybridization is the largest r, sp2 is the largest esd. Mix-n-match.
    covalent_radii[:C] = Dict(:radius_Å=> covalent_radii[:C_sp3][:radius_Å], :esd_Å => covalent_radii[:C_sp2][:esd_Å])
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
COVALENT_RADII = _CorderoParams([:atom, :covalent_radius_A, :esd_Å],
"""H,0.31,.05
He,0.28,.01
Li,1.28,.07
Be,0.96,.03
B,0.84,.03
C_sp3,0.76,.01
C_sp2,0.73,.02
C_sp,0.69,.01
N,0.71,.01
O,0.66,.02
F,0.57,.03
Ne,0.58,.01
Na,1.66,.09
Mg,1.41,.07
Al,1.21,.04
Si,1.11,.02
P,1.07,.03
S,1.05,.03
Cl,1.02,.04
Ar,1.06,.010
K,2.03,.012
Ca,1.76,.010
Sc,1.7,.07
Ti,1.6,.08
V,1.53,.08
Cr,1.39,.05
Mn_lo,1.39,.05
Mn_hi,1.61,.08
Fe_lo,1.32,.03
Fe_hi,1.52,.06
Co_lo,1.26,.03
Co_hi,1.5,.07
Ni,1.24,.04
Cu,1.32,.04
Zn,1.22,.04
Ga,1.22,.03
Ge,1.2,.04
As,1.19,.04
Se,1.2,.04
Br,1.2,.03
Kr,1.16,.04
Rb,2.2,.09
Sr,1.95,.010
Y,1.9,.07
Zr,1.75,.07
Nb,1.64,.06
Mo,1.54,.05
Tc,1.47,.07
Ru,1.46,.07
Rh,1.42,.07
Pd,1.39,.06
Ag,1.45,.05
Cd,1.44,.09
In,1.42,.05
Sn,1.39,.04
Sb,1.39,.05
Te,1.38,.04
I,1.39,.03
Xe,1.4,.09
Cs,2.44,.011
Ba,2.15,.011
La,2.07,.08
Ce,2.04,.09
Pr,2.03,.07
Nd,2.01,.06
Pm,1.99,.01
Sm,1.98,.08
Eu,1.98,.06
Gd,1.96,.06
Tb,1.94,.05
Dy,1.92,.07
Ho,1.92,.07
Er,1.89,.06
Tm,1.9,.010
Yb,1.87,.08
Lu,1.87,.08
Hf,1.75,.010
Ta,1.7,.08
W,1.62,.07
Re,1.51,.07
Os,1.44,.04
Ir,1.41,.06
Pt,1.36,.05
Au,1.36,.06
Hg,1.32,.05
Tl,1.45,.07
Pb,1.46,.05
Bi,1.48,.04
Po,1.4,.04
At,1.5,.01
Rn,1.5,.01
Fr,2.6,.01
Ra,2.21,.02
Ac,2.15,.01
Th,2.06,.06
Pa,2,.01
U,1.96,.07
Np,1.9,.01
Pu,1.87,.01
Am,1.8,.06
Cm,1.69,.03
""")
