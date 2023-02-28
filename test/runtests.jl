using Test
using SafeTestsets: @safetestset

MODULES = [
    # Input modules
    "materials.jl", "cross_sections.jl", "elements.jl", "boundary_conditions.jl", "meshes.jl",
    # Analysis modules
    "structural_model.jl", "structural_solvers.jl", "static_analysis.jl"]

EXAMPLES = [
    joinpath("..", "examples", "von_misses_truss", "von_misses_truss.jl"),
    joinpath("..", "examples", "uniaxial_extension", "uniaxial_extension.jl")
]

function test(files::Vector)
    @testset "ONSAS.jl" begin
        foreach(test, files)
    end
end

function test(file::String)
    @info "Testing $file..."
    path = joinpath(@__DIR__, file)
    @eval @time @safetestset $file begin
        include($path)
    end
end

test([MODULES; EXAMPLES]);
