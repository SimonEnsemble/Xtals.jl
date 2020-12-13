# makes a dataframe from stringified atomicmasses.csv
function _AtomicMasses(headers::Array{Symbol}, data::String)
    species = Symbol[]
    mass = Float64[]
    for line ∈ split(data,"\n")
        fields = split(line, ",")
        if length(fields) == 2
            push!(species, Symbol(fields[1]))
            push!(mass, parse(Float64, fields[2]))
        end
    end
	atomic_masses = Dict{Symbol, Float64}()
    for row in eachrow(DataFrame([species, mass], headers))
		atomic_masses[Symbol(row[:atom])] = row[:mass]
    end
    return atomic_masses
end


# return the global atomic masses dataframe
function get_atomic_masses()
	return ATOMIC_MASSES
end


# Global atomic masses dataframe
ATOMIC_MASSES = _AtomicMasses([:atom, :mass],
"""VOID,0
re,0
He,4.0026
N_in_N2,14.0067
CH2,14.025
CH3,15.035
CH4,16.04
dummy,28.05
Rn,222
Ar,39.948
Ag,107.8682
Al,26.981539
As,74.9216
Au,196.96657
B,10.811
Ba,137.327
Be,9.012182
Bi,208.9804
Br,79.904
C,12.0107
C_b,12.0107
C_tol,12.0107
C_ac,12.0107
C_RCOO,12.0107
C_sp2,12.0107
C_sp3,12.0107
C_CO2,12.0107
Ca,40.078
Cd,112.411
Ce,140.116
Cl,35.453
Co,58.933195
Cr,51.996
Cs,132.9054
Cu,63.546
Dy,162.5
Er,167.259
Eu,151.964
F,18.998403
Fe,55.845
Ga,69.723
Gd,157.25
Ge,72.64
H,1.00794
H_H2S,1.00794
H_b,1.00794
Hf,178.49
Hg,200.59
Ho,164.93032
I,126.90447
In,114.818
Ir,192.217
K,39.0983
Kr,83.798
La,138.9055
Li,6.941
Lu,174.967
Mg,24.305
Mn,54.938
Mo,95.94
N,14.0067
Na,22.989769
Nb,92.90638
Nd,144.242
Ni,58.6934
Np,237
O,15.9994
O_RCOO,15.9994
O_CO2,15.9994
O_zeo,15.9994
P,30.973762
Pb,207.2
Pd,106.42
Pr,140.90765
Pt,195.084
Rb,85.4678
Re,186.207
Rh,102.9055
Ru,101.07
S,32.065
S_H2S,32.065
Sb,121.76
Sc,44.955912
Se,78.96
Si,28.0855
Si_zeo,28.0855
Sm,150.36
Sn,118.71
Sr,87.62
Tb,158.92535
Te,127.6
Th,232.03806
Ti,47.867
Tm,168.93421
V,50.941
Y,88.90585
Yb,173.04
Xe,131.293
Xe_,131.293
Xe__,131.293
Zn,65.38
Zr,91.224
ig,1.0
X,1
""")
