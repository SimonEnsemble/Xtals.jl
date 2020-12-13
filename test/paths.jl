module Path_Test

using Xtals
using Test

@testset "Path Tests" begin
    set_path_to_data(Xtals.PATH_TO_DATA, print=true)
    set_path_to_crystals(joinpath(pwd(), "other_data", "other_crystals"))
    print_file_paths()
    Crystal("other_SBMOF-1.cif")
    @test true # if made it this far :)
end
end
