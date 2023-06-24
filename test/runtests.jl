using Test
using SafeTestsets: @safetestset

MODULES = ["boundary_conditions/boundary_conditions.jl",
           "cross_sections/cross_sections.jl",
           "entities/nodes.jl",
           "entities/tetrahedrons.jl",
           "entities/triangular_faces.jl",
           "entities/trusses.jl",
           "interfaces/gmsh.jl",
           "interfaces/vtk.jl",
           "materials/materials.jl",
           "meshes/handlers.jl",
           "meshes/interpolators.jl",
           "meshes/meshes.jl",
           "meshes/searches.jl",
           "structural_analyses/static_analyses.jl",
           "structural_model/structures.jl",
           "structural_solvers/structural_solvers.jl",
           "utils.jl"]

EXAMPLES_FOLDER = joinpath("..", "examples")

EXAMPLE_NAMES = ["von_misses_truss",
                 "linear_extension",
                 "uniaxial_extension",
                 "uniaxial_compression",
                 "cylinder_internal_pressure"]

EXAMPLES = [joinpath(EXAMPLES_FOLDER, name, name * ".jl") for name in EXAMPLE_NAMES]

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
