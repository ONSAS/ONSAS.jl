"""
Module to handle Dirichlet boundary conditions.
"""
module DirichletBoundaryConditions

using ..Utils, ..BoundaryConditions

export Dirichlet

""" Generalized displacement boundary condition struct.
### Fields:
- `dofs`    -- Vectors of symbols the where the boundary condition is subscripted. 
- `values`  -- Values imposed function. 
- `name`    -- Boundary condition label.
"""
Base.@kwdef struct Dirichlet <: AbstractDirichletBoundaryCondition
    dofs::Vector{Symbol} = [:u]
    values::Function
    name::Label = NO_LABEL
end

end
