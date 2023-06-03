
"""
Fixed degrees of freedom boundary conditions.
"""
module FixedDofBoundaryConditions

using Reexport

using ..BoundaryConditions
using ..Nodes
using ..Entities
using ..Utils

@reexport import ..BoundaryConditions: apply

export FixedDof, components

"""
Fixed displacement boundary condition.

This is a particular instance of the struct `DisplacementBoundaryCondition`
considering null dof value at an specific component of the dof displacements.
"""
Base.@kwdef struct FixedDof <: AbstractDirichletBoundaryCondition
    "Symbols where the where the boundary condition is subscribed."
    dofs::Vector{Field} = [:u]
    "Vectors of integer indicating the fixed degree of freedom component."
    components::Vector{Dof}
    "Label of the boundary condition."
    name::Label = NO_LABEL
end

"Return the fixed components of the `Dof`s defined in the boundary condition `bc`."
components(bc::FixedDof) = bc.components

"Return fixed `Dof`s of an `AbstractNode` imposed in the `FixedDof` `fbc`."
function apply(fbc::FixedDof, n::AbstractNode)
    fbc_dofs_symbols = dofs(fbc)
    dofs_to_delete = Dof[]
    for dof_symbol in fbc_dofs_symbols
        push!(dofs_to_delete, getindex(dofs(n), dof_symbol)[components(fbc)]...)
    end
    dofs_to_delete
end

"Return fixed `Dof`s of an `AbstractFace` or `AbstractElement` imposed in the `FixedDof` `fbc`."
function apply(fbc::FixedDof, e::E) where {E<:Union{AbstractFace,AbstractElement}}
    dofs_to_delete = Dof[]
    [push!(dofs_to_delete, apply(fbc, n)...) for n in nodes(e)]
    dofs_to_delete
end

end
