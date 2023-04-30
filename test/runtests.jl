using Test
using SafeTestsets: @safetestset

MODULES = ["materials.jl",
           "cross_sections.jl",
           "elements.jl",
           "boundary_conditions.jl",
           "meshes.jl",
           "interfaces/gmsh.jl",
           "interfaces/vtk.jl",
           "structural_model.jl",
           "structural_solvers.jl",
           "static_analyses.jl"]

EXAMPLES = [joinpath("..", "examples", "von_misses_truss", "von_misses_truss.jl"),
            joinpath("..", "examples", "linear_extension", "linear_extension.jl"),
            joinpath("..", "examples", "uniaxial_extension", "uniaxial_extension.jl"),
            joinpath("..", "examples", "uniaxial_compression", "uniaxial_compression.jl")]

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

test([MODULES, EXAMPLES]);
