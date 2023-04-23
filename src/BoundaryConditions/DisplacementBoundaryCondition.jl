using ..BoundaryConditions: AbstractDisplacementBoundaryCondition
using ..Elements: AbstractNode, AbstractFace, AbstractElement

export DisplacementBoundaryCondition

""" Generalized displacement boundary condition struct.
### Fields:
- `dofs`    -- Vectors of symbols the where the boundary condition is subscripted. 
- `values`  -- Values imposed function. 
- `name`    -- Boundary condition label.
"""
Base.@kwdef struct DisplacementBoundaryCondition <: AbstractDisplacementBoundaryCondition
    dofs::Vector{Symbol}
    values::Function
    name::Symbol = :no_labelled_bc
end
