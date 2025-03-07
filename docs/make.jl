using Documenter
using Pkg: Pkg

# Fix for https://github.com/trixi-framework/Trixi.jl/issues/668
# to allow building the docs locally
if (get(ENV, "CI", nothing) != "true") &&
   (get(ENV, "JULIA_DOC_DEFAULT_ENVIRONMENT", nothing) != "true")
    push!(LOAD_PATH, dirname(@__DIR__))
end

using RKOpt

# Define module-wide setups such that the respective modules
# are available in doctests
DocMeta.setdocmeta!(RKOpt,
                    :DocTestSetup,
                    :(using RKOpt);
                    recursive = true)

# Copy some files from the top level directory to the docs and modify them
# as necessary
open(joinpath(@__DIR__, "src", "license.md"), "w") do io
    # Point to source license file
    println(io, """
    ```@meta
    EditURL = "https://github.com/NumericalMathematics/RKOpt.jl/blob/main/LICENSE"
    ```
    """)
    # Write the modified contents
    println(io, "# License")
    println(io, "")
    for line in eachline(joinpath(dirname(@__DIR__), "LICENSE.md"))
        line = replace(line, "[LICENSE.md](LICENSE.md)" => "[License](@ref)")
        println(io, "> ", line)
    end
end

open(joinpath(@__DIR__, "src", "code_of_conduct.md"), "w") do io
    # Point to source license file
    println(io,
            """
            ```@meta
            EditURL = "https://github.com/NumericalMathematics/RKOpt.jl/blob/main/CODE_OF_CONDUCT.md"
            ```
            """)
    # Write the modified contents
    println(io, "# [Code of Conduct](@id code-of-conduct)")
    println(io, "")
    for line in eachline(joinpath(dirname(@__DIR__), "CODE_OF_CONDUCT.md"))
        line = replace(line, "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref)")
        println(io, "> ", line)
    end
end

open(joinpath(@__DIR__, "src", "contributing.md"), "w") do io
    # Point to source license file
    println(io, """
    ```@meta
    EditURL = "https://github.com/NumericalMathematics/RKOpt.jl/blob/main/CONTRIBUTING.md"
    ```
    """)
    # Write the modified contents
    println(io, "# Contributing")
    println(io, "")
    for line in eachline(joinpath(dirname(@__DIR__), "CONTRIBUTING.md"))
        line = replace(line, "[LICENSE.md](LICENSE.md)" => "[License](@ref)")
        println(io, "> ", line)
    end
end

# Make documentation
makedocs(modules = [RKOpt],
         sitename = "RKOpt.jl",
         format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true",
                                  canonical = "https://NumericalMathematics.github.io/RKOpt.jl/stable"),
         # Explicitly specify documentation structure
         pages = [
             "Home" => "index.md",
             "API reference" => "api_reference.md",
             "Contributing" => "contributing.md",
             "Code of conduct" => "code_of_conduct.md",
             "License" => "license.md"
         ])

deploydocs(repo = "github.com/NumericalMathematics/RKOpt.jl",
           devbranch = "main",
           push_preview = true)
