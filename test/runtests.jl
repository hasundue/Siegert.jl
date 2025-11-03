# Top-level test runner; include split test files

test_files = filter(f -> startswith(f, "test_"), readdir(@__DIR__))

if isempty(ARGS)
    for test_file in test_files
        include(test_file)
    end
else
    for arg in ARGS
        if arg in test_files
            include(arg)
        elseif startswith(arg, "test_") && endswith(arg, ".jl")
            @warn "Test file not found: $arg"
        else
            matched = filter(f -> occursin(arg, f), test_files)
            if isempty(matched)
                @warn "No test file matching: $arg"
            else
                for test_file in matched
                    include(test_file)
                end
            end
        end
    end
end
