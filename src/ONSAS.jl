module ONSAS

using Reexport: @reexport

FILES = ["Utils.jl",
         "Materials/Materials.jl",
         "CrossSections/CrossSections.jl",
         "Elements/Elements.jl",
         "BoundaryConditions/BoundaryConditions.jl",
         "Meshes/Meshes.jl",
         "Meshes/Handlers.jl",
         "Interfaces/Gmsh.jl",
         "Interfaces/VTK.jl",
         "StructuralModel/StructuralModel.jl",
         "StructuralSolvers/StructuralSolvers.jl",
         "StructuralAnalyses/StructuralAnalyses.jl"]

foreach(FILES) do m
    include(m)
end

# Utility methods.
@reexport using .Utils

# Physical models.
@reexport using .Materials
@reexport using .CrossSections
@reexport using .Elements
@reexport using .BoundaryConditions

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

end # module
