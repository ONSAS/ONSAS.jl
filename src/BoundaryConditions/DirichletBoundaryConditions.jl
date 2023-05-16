"""
Module to handle Dirichlet boundary conditions.
"""
module DirichletBoundaryConditions

using ..Utils
using ..BoundaryConditions

export Dirichlet

"""
Generalized displacement boundary condition. 
"""
Base.@kwdef struct Dirichlet <: AbstractDirichletBoundaryCondition
    "Vector of fields for which this boundary condition applies to."
    dofs::Vector{Field} = [:u]
    "Function defining values which are imposed by this boundary condition."
    values::Function
    "Boundary condition label."
    name::Label = NO_LABEL
end

end
