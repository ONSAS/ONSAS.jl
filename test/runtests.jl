using SafeTestsets: @safetestset

@safetestset "ONSAS.Mesh" begin
    include("meshes.jl")
end

@safetestset "ONSAS.Materials" begin
    include("materials.jl")
end

@safetestset "ONSAS.BoundaryConditions" begin
    include("boundary_conditions.jl")
end

@safetestset "ONSAS.Elements" begin
    include("elements.jl")
end

@safetestset "ONSAS.StructuralModel" begin
    include("structural_model.jl")
end
