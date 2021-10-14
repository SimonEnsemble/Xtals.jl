using Test, LightGraphs, MetaGraphs, Documenter, FIGlet

if !isdir("temp")
    mkdir("temp")
end

FIGlet.render("Xtals.jl", FIGlet.availablefonts()[5])

using Xtals

# run doctests
doctest(Xtals)

@info "Tests complete!"
