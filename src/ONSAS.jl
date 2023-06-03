module ONSAS

using Reexport: @reexport

FILES = ["Utils.jl",
         # Materials
         "Materials/Materials.jl",
         "Materials/LinearElasticMaterials.jl",
         "Materials/HyperElasticMaterials.jl",
         "Materials/IsotropicLinearElasticMaterial.jl",
         "Materials/SvkMaterial.jl",
         "Materials/NeoHookeanMaterial.jl",
         "Materials/HyperElasticMaterial.jl",
         # Cross-sections
         "CrossSections/CrossSections.jl",
         "CrossSections/Circles.jl",
         "CrossSections/Rectangles.jl",
         "CrossSections/Squares.jl",
         "CrossSections/GenericCrossSections.jl",
         # Elements
         "Elements/Elements.jl",
         # Boundary conditions
         "BoundaryConditions/BoundaryConditions.jl",
         "BoundaryConditions/FixedDofBoundaryConditions.jl",
         "BoundaryConditions/DirichletBoundaryConditions.jl",
         "BoundaryConditions/GlobalLoadBoundaryConditions.jl",
         "BoundaryConditions/LocalLoadBoundaryConditions.jl",
         # Meshes
         "Meshes/Meshes.jl",
         "Meshes/Interpolators.jl",
         "Meshes/Handlers.jl",
         # Interfaces
         "Interfaces/Gmsh.jl",
         "Interfaces/VTK.jl",
         #Structural Model
         "StructuralModel/StructuralModel.jl",
         # Structural Solvers
         "StructuralSolvers/StructuralSolvers.jl",
         "StructuralSolvers/Solvers.jl",
         "StructuralAnalyses/StructuralAnalyses.jl",
         "StructuralAnalyses/StaticAnalyses.jl"]

foreach(FILES) do m
    include(m)
end

# Utility methods.
@reexport using .Utils

# Physical models.
## Materials
@reexport using .Materials
@reexport using .LinearElasticMaterials
@reexport using .IsotropicLinearElasticMaterial
@reexport using .HyperElasticMaterials
@reexport using .SvkMaterial
@reexport using .NeoHookeanMaterial
@reexport using .HyperElasticMaterial

## Cross-sections
@reexport using .CrossSections
@reexport using .Circles
@reexport using .Rectangles
@reexport using .Squares
@reexport using .GenericCrossSections

## Entities
@reexport using .Elements

## Boundary conditions
@reexport using .BoundaryConditions
@reexport using .FixedDofBoundaryConditions
@reexport using .DirichletBoundaryConditions
@reexport using .LocalLoadBoundaryConditions
@reexport using .GlobalLoadBoundaryConditions

# Geometric entities and interpolation.
@reexport using .Meshes
@reexport using .Interpolators
@reexport using .Handlers

# Interfaces with external programs.
@reexport using .Gmsh
@reexport using .VTK

# Structural models.
@reexport using .StructuralModel

# Finite element solvers.
@reexport using .StructuralSolvers
@reexport using .Solvers

# Structural analyses.
@reexport using .StructuralAnalyses
@reexport using .StaticAnalyses

end
