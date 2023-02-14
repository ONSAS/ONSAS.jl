module ONSAS

# load packages
include("init.jl")


# Input modules
include("interface/Materials.jl")
@reexport using .Materials

include("interface/BoundaryConditions.jl")
@reexport using .BoundaryConditions

include("interface/CrossSections.jl")
@reexport using .CrossSections

include("interface/Elements.jl")
@reexport using .Elements

include("interface/Meshes.jl")
@reexport using .Meshes

# Structural model 
include("interface/StructuralModel.jl")
@reexport using .StructuralModel

# # Solvers
# include("interface/StructuralSolvers.jl")
# @reexport using .StructuralSolvers

# # Analysis
# include("interface/StructuralAnalyses.jl")
# @reexport using .StructuralAnalyses


end # module
