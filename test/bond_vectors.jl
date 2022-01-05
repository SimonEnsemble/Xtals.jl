@testset "bond vectors" begin
    box = replicate(unit_cube(), (10, 10, 10))
    xtal = Crystal(
        "vector test 1",
        box,
        Atoms(
            [:H, :H, :H],
            Frac([
                0.000 0.000 0.075;
                0.075 0.000 0.000;
                0.000 0.000 0.000
            ])
        ),
        Charges{Frac}(0)
    )

    # check non-fatal error for calculating vectors w/o bonds
    @test isnothing(Xtals.calculate_bond_vectors!(xtal))

    infer_bonds!(xtal, true, calculate_vectors=true)

    # check fatal error for attempting to re-calculate vectors
    @test_throws ErrorException calculate_bond_vectors!(xtal)

    # test clearing and re-calculating vectors
    clear_vectors!(xtal)

    @test !get_prop(xtal.bonds, :has_vectors)

    calculate_bond_vectors!(xtal)

    @test get_prop(xtal.bonds, :has_vectors)

    # check vectors
    @test isapprox(get_bond_vector(xtal, 1, 2), [0.00; -0.75;  0.00])

    @test isapprox(get_bond_vector(xtal, 2, 3), [0.75;  0.00;  0.00])

    # check reversal of direction
    @test isapprox(get_bond_vector(xtal, 2, 1), [0.00;  0.75;  0.00])

    # check bond angle
    @test isapprox(bond_angle(xtal, 1, 2, 3), deg2rad(90))

    # check invalid bond angle
    @test_throws AssertionError bond_angle(xtal, 2, 3, 1)

    # check bond distance calculations
    @test isapprox(bond_distance(xtal, 1, 2), norm(get_prop(xtal.bonds, 1, 2, :vector)))

    xtal = Crystal(
        "vector test 2",
        box,
        Atoms(
            [:H, :H, :H],
            Frac([
                1.000 0.000 0.075;
                0.075 0.000 0.000;
                0.000 0.000 0.000
            ])
        ),
        Charges{Frac}(0)
    )
    infer_bonds!(xtal, true, calculate_vectors=true)

    @test isapprox(get_bond_vector(xtal, 1, 2), [0.00; -0.75;  0.00])

    @test isapprox(get_bond_vector(xtal, 2, 1), [0.00;  0.75;  0.00])
    
    @test isapprox(get_bond_vector(xtal, 2, 3), [0.75;  0.00;  0.00])
    
    @test isapprox(bond_angle(xtal, 1, 2, 3), deg2rad(90)) # 90° angles work to numerical precision

    xtal = Crystal(
        "vector test 3",
        box,
        Atoms(
            [:H, :H, :H],
            Frac([
                0.000 0.075 0.150;
                0.000 0.000 0.000;
                0.000 0.000 0.000
            ])
        ),
        Charges{Frac}(0)
    )
    infer_bonds!(xtal, true, calculate_vectors=true)
    @test isapprox(get_bond_vector(xtal, 1, 2), [0.75; 0.00; 0.00])
    
    @test isapprox(get_bond_vector(xtal, 2, 3), [0.75; 0.00; 0.00])
    
    @test isapprox(bond_angle(xtal, 1, 2, 3), π, rtol=1e-4) # 180° angles work to within 0.01%

    methane = Frac(read_xyz(joinpath(rc[:paths][:crystals], "methane.xyz")), box)
    xtal = Crystal(
        "vector test 4",
        box,
        methane,
        Charges{Frac}(0)
    )
    infer_bonds!(xtal, true, calculate_vectors=true)
    
    # test that ∠ijk = ∠kji and all angles are all the same
    α = bond_angle(xtal, 2, 1, 3)
    results = Bool[]
    for i ∈ 2:5
        for j ∈ 2:5
            if i == j
                push!(results, isapprox(bond_angle(xtal, i, 1, j), 0)) # 0° angles work to numerical precision
            else
                push!(results, isapprox(bond_angle(xtal, i, 1, j), α, rtol=1e-5)) # identical angles work to within 0.001%
            end
        end        
    end
    
    @test all(results)
    
    # test that bond angles match expected for tetrahedral carbon
    @test isapprox(bond_angle(xtal, 2, 1, 3), deg2rad(109.5), rtol=1e-3) # 109.5° angles work to within 0.1%
    
    # test that result is translation invariant
    results = Vector{Bool}(undef, 10000)
    for i ∈ 1:length(results)
        translate_by!(xtal.atoms.coords, Frac(rand(3, 1)))
        results[i] = isapprox(bond_angle(xtal, 2, 1, 3), α, rtol=1e-5)
    end
    
    @test all(results)
end