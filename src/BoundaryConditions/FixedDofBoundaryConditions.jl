
"""
Fixed degrees of freedom boundary conditions.
"""
module FixedDofBoundaryConditions

using Reexport
using ..BoundaryConditions, ..Elements, ..Utils

@reexport import ..BoundaryConditions: apply

export components

""" Fixed displacement boundary condition struct:
This is a particular instance of the struct `DisplacementBoundaryCondition`
    considering null dof value at an specific component of the dof displacements
"""
Base.@kwdef struct FixedDof <: AbstractDirichletBoundaryCondition
    "Symbols where the where the boundary condition is subscribed."
    dofs::Vector{Symbol} = [:u]
    "Vectors of integer indicating the degree of freedom component fixed."
    components::Vector{Int}
    "Label of the boundary condition."
    name::Symbol = NO_LABEL
end

"Constructor of `FixedDof` boundary condition given a `Vector` of `Dof`s symbols and the fixed `components`."
function FixedDof(dofs::Vector{Symbol}, components::Vector{Int},
                  name::String=string(NO_LABEL))
    FixedDof(dofs, components, Symbol(name))
end

"Return the fixed components of the `Dof`s defined in the boundary condition `bc`."
components(bc::FixedDof) = bc.components

"Return fixed `Dof`s of an `AbstractNode` imposed in the `FixedDof` `fbc`."
function apply(fbc::FixedDof, n::AbstractNode)
    fbc_dofs_symbols = dofs(fbc)
    dofs_to_delete = Vector{Dof}(undef, length(fbc_dofs_symbols))
    for dof_symbol in fbc_dofs_symbols
        push!(dofs_to_delete, getindex(dofs(n), dof_symbol)[components(fbc)]...)
    end
    dofs_to_delete
end

"Return fixed `Dof`s of an `AbstractFace` or `AbstractElement` imposed in the `FixedDof` `fbc`."
function apply(fbc::FixedDof, e::E) where {E<:Union{AbstractFace,AbstractElement}}
    dofs_to_delete = Dof[]
    [push!(dofs_to_delete, _apply(fbc, n)...) for n in nodes(e)]
    dofs_to_delete
end

end
