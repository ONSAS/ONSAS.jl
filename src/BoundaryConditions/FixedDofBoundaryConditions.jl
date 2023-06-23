
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
Fixed boundary condition.

Considers null dof values at specific component(s) of the given field.
"""
struct FixedDof <: AbstractDirichletBoundaryCondition
    "Field where the boundary condition applies to."
    field::Field
    "Components of the field which are fixed."
    components::Vector{Int64}
    "Boundary condition label."
    name::Label
    function FixedDof(field::Field, components::Vector{Int64}, name::Label=NO_LABEL)
        new(field, components, name)
    end
end

"Return the fixed components of the `Dof`s defined in the boundary condition `bc`."
components(bc::FixedDof) = bc.components

"Return fixed `Dof`s of an `AbstractNode` imposed in the `FixedDof` `fbc`."
function apply(bc::FixedDof, n::AbstractNode)
    # TODO Rename method to fixed_dofs ?
    dofs(n, bc.field)[bc.components]
end

"Return fixed `Dof`s of an `AbstractFace` or `AbstractElement` imposed in the `FixedDof` `fbc`."
function apply(bc::FixedDof, e::E) where {E<:Union{AbstractFace,AbstractElement}}
    # TODO Rename method to fixed_dofs ?
    reduce(vcat, apply(bc, n) for n in nodes(e))
end

end
