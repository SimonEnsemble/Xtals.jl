module Path_Test

using Xtals
using Test

@testset "Path Tests" begin
    set_global(:path_to_crystals => joinpath(pwd(), "other_data", "other_crystals"))
    print_file_paths()
    Crystal("other_SBMOF-1.cif")
    @test true # if made it this far :)
end
end
