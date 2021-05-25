"""
    radius_dict = parse_covalent_radii(data)

Create a dictionary with the covalent radius for each element, using data from Cordero (10.1039/B801115J) and Manz (10.1039/c9ra07327b)

# Returns
- `radius_dict::Dict{Symbol,Float64}`: A dictionary with elements as keys and covalent radii (Å) as the values.
"""
function parse_covalent_radii(data::String)::Dict{Symbol,Float64}
    radius_dict = Dict()
    for line ∈ split(data,"\n")
        fields = split(line, ",")
        if length(fields) == 2
            radius_dict[Symbol(fields[1])] = parse(Float64, fields[2])
        end
    end
    return radius_dict
end


# Cordero bond radii, DOI: 10.1039/B801115J Table 2.
# Overridden w/ larger values from DOI 10.1039/c9ra07327b Table 5
global DEFAULT_COVALENT_RADII = parse_covalent_radii(
    """
    H,0.38
    He,0.28
    Li,1.28
    Be,0.96
    B,1.01
    C,0.88
    N,0.86
    O,0.89
    F,0.82
    Ne,0.58
    Na,1.66
    Mg,1.41
    Al,1.53
    Si,1.38
    P,1.28
    S,1.20
    Cl,1.17
    Ar,1.06
    K,2.03
    Ca,1.76
    Sc,1.7
    Ti,1.65
    V,1.53
    Cr,1.39
    Mn,1.61
    Fe,1.52
    Co,1.5
    Ni,1.24
    Cu,1.32
    Zn,1.22
    Ga,1.22
    Ge,1.2
    As,1.19
    Se,1.2
    Br,1.2
    Kr,1.16
    Rb,2.2
    Sr,1.95
    Y,1.9
    Zr,1.75
    Nb,1.64
    Mo,1.54
    Tc,1.47
    Ru,1.46
    Rh,1.42
    Pd,1.68
    Ag,1.56
    Cd,1.56
    In,1.53
    Sn,1.64
    Sb,1.64
    Te,1.65
    I,1.58
    Xe,1.4
    Cs,2.44
    Ba,2.15
    La,2.07
    Ce,2.04
    Pr,2.03
    Nd,2.01
    Pm,1.99
    Sm,1.98
    Eu,1.98
    Gd,1.96
    Tb,1.94
    Dy,1.92
    Ho,1.92
    Er,1.89
    Tm,1.9
    Yb,1.87
    Lu,1.87
    Hf,1.75
    Ta,1.7
    W,1.62
    Re,1.51
    Os,1.44
    Ir,1.50
    Pt,1.66
    Au,1.68
    Hg,1.88
    Tl,1.45
    Pb,1.72
    Bi,1.72
    Po,1.4
    At,1.5
    Rn,1.5
    Fr,2.6
    Ra,2.21
    Ac,2.15
    Th,2.06
    Pa,2.0
    U,1.96
    Np,1.9
    Pu,1.87
    Am,1.8
    Cm,1.69
    """
)
