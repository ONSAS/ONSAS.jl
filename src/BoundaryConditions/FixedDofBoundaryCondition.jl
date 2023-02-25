using ..BoundaryConditions: AbstractDisplacementBoundaryCondition
using ..Elements: Dof, AbstractNode, AbstractFace, AbstractElement, nodes

import ..BoundaryConditions: _apply

export components

""" Fixed displacement boundary condition struct:
This is a particular instance of the struct `DisplacementBoundaryCondition`
    considering null dof value at an specific component of the dof displacements
### Fields:
- `dofs`             -- Vectors of symbols the where the boundary condition is subscribed.
- `components` -- Vectors of integer indicating the degree of freedom component fixed.
- `name`             -- Boundary condition label.
"""
Base.@kwdef struct FixedDofBoundaryCondition <: AbstractDisplacementBoundaryCondition
    dofs::Vector{Symbol}
    components::Vector{Int}
    name::Symbol = :no_labelled_bc
end

FixedDofBoundaryCondition(dofs::Vector{Symbol}, components::Vector{Int}, name::String="no_labelled_bc") =
    FixedDofBoundaryCondition(dofs, components, Symbol(name))

"Returns the fixed components of the `Dof`s defined in the boundary condition `bc`."
components(bc::FixedDofBoundaryCondition) = bc.components


"Returns fixed `Dof`s of an `AbstractNode` imposed in the `FixedDofBoundaryCondition` `fbc`."
function _apply(fbc::FixedDofBoundaryCondition, n::AbstractNode)
    fbc_dofs_symbols = dofs(fbc)
    dofs_to_delete = Dof[]
    for dof_symbol in fbc_dofs_symbols
        push!(dofs_to_delete, getindex(dofs(n), dof_symbol)[components(fbc)]...)
    end
    dofs_to_delete
end

"Returns fixed `Dof`s of an `AbstractFace` or `AbstractElement` imposed in the `FixedDofBoundaryCondition` `fbc`."
function _apply(fbc::FixedDofBoundaryCondition, e::E) where {E<:Union{AbstractFace,AbstractElement}}
    dofs_to_delete = Dof[]
    [push!(dofs_to_delete, _apply(fbc, n)...) for n in nodes(e)]
    dofs_to_delete
end
