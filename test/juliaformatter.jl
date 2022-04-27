using DeepDiffs, JuliaFormatter

#check_dirs = ["src", "test"]
check_dirs = ["."]

function formatting_diff(filename::String)
    @assert contains(filename, ".jl")
    @assert isfile(filename) "Cannot access $filename"
    # read file as string
    original = String(Char.(read(filename)))
    # format file string
    formatted = format_text(original)
    # compare
    return deepdiff(original, formatted)
end


function check_file(file::String)::Bool
    diff = formatting_diff(file)
    if before(diff) ≠ after(diff)
        char_idx = removed(diff)
        lines = []
        for i in char_idx
            try
                push!(lines, length(findall("\n", before(diff)[1:i])) + 1)
            catch exception
                @warn file exception
            end
        end
        lines = chop(reduce(*, ["$x, " for x in unique(lines)]), tail=2)
        @error "Possible formatting issue" file lines
        return false
    end
    return true
end


function check_all(check_dirs::Vector{String})::Bool
    no_fails = true
    for dir in check_dirs
        for file in [joinpath(dir, file) for file in readdir(dir) if contains(file, ".jl")]
            if !check_file(file)
                no_fails = false
            end
        end
    end
    return no_fails
end

@assert check_all(check_dirs) "Possible formatting issues detected."
