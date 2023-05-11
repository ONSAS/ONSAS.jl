module ONSAS

using Reexport: @reexport

FILES = ["Utils.jl",
         "Materials/Materials.jl",
         "CrossSections/CrossSections.jl",
         "Elements/Elements.jl",
         # Boundary conditions
         "BoundaryConditions/BoundaryConditions.jl",
         "BoundaryConditions/FixedDofBoundaryConditions.jl",
         "BoundaryConditions/GlobalLoadBoundaryConditions.jl",
         "BoundaryConditions/LocalLoadBoundaryConditions.jl",
         # Meshes
         "Meshes/Meshes.jl",
         "Meshes/Handlers.jl",
         # Interfaces
         "Interfaces/Gmsh.jl",
         "Interfaces/VTK.jl",
         #Structural Model
         "StructuralModel/StructuralModel.jl",
         # Structural Solverss
         "StructuralSolvers/StructuralSolvers.jl",
         "StructuralAnalyses/StructuralAnalyses.jl",
         "StructuralAnalyses/StaticAnalyses.jl"]

foreach(FILES) do m
    include(m)
end

# Utility methods.
@reexport using .Utils

# Physical models.
# Materials
@reexport using .Materials
# Cross-sections
@reexport using .CrossSections
# Entities
@reexport using .Elements
# Boundary conditions
@reexport using .BoundaryConditions
@reexport using .FixedDofBoundaryConditions
@reexport using .DirichletBoundaryConditions
@reexport using .NeumannBoundaryConditions
@reexport using .LocalLoadBoundaryConditions
@reexport using .GlobalLoadBoundaryConditions

#=
# GEometric entities and interpolation.
@reexport using .Meshes
@reexport using .Handlers

# Interfaces with external programs.
@reexport using .Gmsh
@reexport using .VTK

# Structural models.
@reexport using .StructuralModel

# Finite element solvers.
@reexport using .StructuralSolvers

# Structural analyses.
@reexport using .StructuralAnalyses
@reexport using .StaticAnalyses
=#
end
