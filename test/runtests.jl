#using Logging
#global_logger(ConsoleLogger(stdout, Logging.Info))
using Test, Revise
using Xtals

function runtest(testfile::String)
    @info "Testing $(testfile)"
    try
        include(testfile)
    catch exception
        @error "Exception in $(testfile)" exception
    end
end

testfiles = ["box.jl",
             "matter.jl",
             "crystal.jl",
             "distance.jl",
             "misc.jl",
             "assert_p1_symmetry.jl",
             "paths.jl" # last, to not interfere with reading test data
             ]

runtest.(testfiles)
