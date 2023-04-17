using Test
using SafeTestsets: @safetestset

MODULES = [
    # Input modules
    "test_Materials.jl", "test_CrossSections.jl", "test_Elements.jl", "test_BoundaryConditions.jl", "test_Meshes.jl",
    # Analysis modules
    "test_StructuralModel.jl", "test_StructuralSolvers.jl", "test_StaticAnalyses.jl"]

EXAMPLES = [
    joinpath("..", "examples", "von_misses_truss", "von_misses_truss.jl"),
    joinpath("..", "examples", "linear_extension", "linear_extension.jl"),
    joinpath("..", "examples", "uniaxial_extension", "uniaxial_extension.jl"),
    joinpath("..", "examples", "uniaxial_compression", "uniaxial_compression.jl"),
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

test([EXAMPLES; MODULES]);
