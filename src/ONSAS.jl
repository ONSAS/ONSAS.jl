module ONSAS

using Reexport: @reexport

# Utils 
include("./Utils.jl")
@reexport using .Utils

# # Input modules
include("Materials/Materials.jl")
@reexport using .Materials

include("CrossSections/CrossSections.jl")
@reexport using .CrossSections

include("Elements/Elements.jl")
@reexport using .Elements

include("BoundaryConditions/BoundaryConditions.jl")
@reexport using .BoundaryConditions

include("Meshes/Meshes.jl")
@reexport using .Meshes

# Structural model 
include("StructuralModel/StructuralModel.jl")
@reexport using .StructuralModel

# # # Solvers
include("StructuralSolvers/StructuralSolvers.jl")
@reexport using .StructuralSolvers

# # # Analysis
# include("interface/StructuralAnalyses.jl")
# @reexport using .StructuralAnalyses


end # module
