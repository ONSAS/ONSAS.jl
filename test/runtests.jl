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

EXAMPLES_FOLDER = joinpath("..", "examples")

EXAMPLE_NAMES = ["von_misses_truss",
                 "linear_extension",
                 "uniaxial_extension",
                 "uniaxial_compression",
                 "cylinder_internal_pressure"]

EXAMPLES = [joinpath(EXAMPLES_FOLDER, EXAMPLE_NAME, EXAMPLE_NAME * ".jl")
            for EXAMPLE_NAME in EXAMPLE_NAMES]

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
