module symmetry_test

using Xtals
using Test

# Test set for making sure simulations meet certain rules
#   i.e. frameworks musts be in P1 symmetry to be used in GCMC or Henry
#   coefficient test
@testset "P1 Symmetry Tests" begin
    non_P1_framework = Crystal("symmetry_test_structure.cif", convert_to_p1=false)
    # Test that an assertion is thrown when trying to replicate a non-P1
    #   structure
    @test_throws ErrorException replicate(non_P1_framework, (2, 2, 2))
    @test_throws ErrorException replicate(non_P1_framework, (1, 1, 1))

end
end
