module ONSAS

using Reexport: @reexport

FILES = ["Utils.jl",
         "Materials/Materials.jl",
         "CrossSections/CrossSections.jl",
         "Elements/Elements.jl",
         "BoundaryConditions/BoundaryConditions.jl",
         "Meshes/Meshes.jl",
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
@reexport using .Meshes

# Structural models.
@reexport using .StructuralModel

# Finite element solvers.
@reexport using .StructuralSolvers

# Structural analyses.
@reexport using .StructuralAnalyses

end # module
