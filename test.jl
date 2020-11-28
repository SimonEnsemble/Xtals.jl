### A Pluto.jl notebook ###
# v0.12.11

using Markdown
using InteractiveUtils

# ╔═╡ d9bc81d0-2f6e-11eb-334a-23b76f7604f8
begin
	push!(LOAD_PATH, joinpath(pwd(), "src"))
	using DataFrames, CSV, Xtals, BenchmarkTools, ProgressMeter, PyPlot
end

# ╔═╡ f2d697ee-2f6e-11eb-3ca4-57f4932912a5
cordero_df = CSV.read("data/cordero_simple.csv", DataFrame)

# ╔═╡ e58d6170-2f71-11eb-28fb-4d4fdc53de7a
cordero_df[!, :x] = Symbol.(cordero_df.x)

# ╔═╡ f96bbed0-2f71-11eb-08ff-0b6af702203b
cordero_df

# ╔═╡ 074eeb60-2f6f-11eb-3395-95cf8e7922e7
md"""
**Rules**
- wildcard: if `:*` is within the smallest diatomic bond distance, it's bonded
- each atom in cordero has one rule for each other atom
- rule (i, j) == rule (j, i); only record one
"""

# ╔═╡ ba3a4740-2f76-11eb-388e-bbc4448f023e
begin
	cordero_df[!, :plus3σ] = 3 .* cordero_df.esd_pm ./ 100 .+ cordero_df.r_Ang
	cordero_df[!, :minus3σ] = cordero_df.r_Ang .- 3 .* cordero_df.esd_pm ./ 100
end

# ╔═╡ cb95d90e-2f70-11eb-36a1-6f4628589b57
begin
	rules = BondingRule[]
	# loop over rows 
	for (i, species_i) ∈ enumerate(eachrow(cordero_df))
		# make every x -> y rule, no commutative duplicates
		for (j, species_j) ∈ enumerate(eachrow(cordero_df))
			if j < i
				continue
			end
			# calculate max and min radii for 3 esd's
			rmin_i = species_i.r_Ang - 3*species_i.esd_pm/100
			rmin_j = species_j.r_Ang - 3*species_j.esd_pm/100
			rmax_i = species_i.r_Ang + 3*species_i.esd_pm/100
			rmax_j = species_j.r_Ang + 3*species_j.esd_pm/100
			rmin = rmin_i + rmin_j
			rmax = rmax_i + rmax_j
			push!(rules, BondingRule(species_i.x, species_j.x, rmin, rmax))
		end
	end
	rules
end

# ╔═╡ 767b5520-2f77-11eb-2c79-c96a59144921
begin
	close()
	figure()
	scatter(1:nrow(cordero_df), cordero_df.r_Ang, label="Cordero")
	scatter(1:nrow(cordero_df), cordero_df.plus3σ, label="Xtals.jl Upper")
	scatter(1:nrow(cordero_df), cordero_df.minus3σ, label="Xtals.jl Lower")
	title("Covalent Radii")
	ylabel("r (Å)   ", rotation="horizontal")
	xlabel("Z")
	legend()
	gcf()
end

# ╔═╡ Cell order:
# ╠═d9bc81d0-2f6e-11eb-334a-23b76f7604f8
# ╠═f2d697ee-2f6e-11eb-3ca4-57f4932912a5
# ╠═e58d6170-2f71-11eb-28fb-4d4fdc53de7a
# ╠═f96bbed0-2f71-11eb-08ff-0b6af702203b
# ╟─074eeb60-2f6f-11eb-3395-95cf8e7922e7
# ╠═ba3a4740-2f76-11eb-388e-bbc4448f023e
# ╠═cb95d90e-2f70-11eb-36a1-6f4628589b57
# ╠═767b5520-2f77-11eb-2c79-c96a59144921
