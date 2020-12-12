TRAVIS = true # set false for better UX in Atom. must be true for CI

testfiles = [
    #"box.jl",
    #"matter.jl",
    "crystal.jl",
    #"distance.jl",
    #"misc.jl",
    #"assert_p1_symmetry.jl",
    "bonds.jl",
    #"paths.jl" # last, to not interfere with reading test data
    ]

using Test
using LightGraphs, MetaGraphs
if !TRAVIS
    using Logging, Revise
    global_logger(ConsoleLogger(stdout, Logging.Info))
end
using Xtals

function runtest(testfile::String)
    @info "Testing $(testfile)"
    if !TRAVIS
        try
            include(testfile)
        catch exception
            @error "Exception in $(testfile)" exception
        end
    else
        include(testfile)
    end
end

if !TRAVIS
    home = dirname(pathof(Xtals))
    set_path_to_data(joinpath(home, "../test/data"),
        relpath_xtals=true, print=true)
    cd(joinpath(home, "../test"))
end

runtest.(testfiles)
