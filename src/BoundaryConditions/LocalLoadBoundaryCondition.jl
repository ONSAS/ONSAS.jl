using ..BoundaryConditions: AbstractLoadBoundaryCondition
using ..Elements: AbstractNode, AbstractFace, AbstractElement, normal_direction

import ..BoundaryConditions: _apply

export LocalLoadBoundaryCondition

""" Load boundary condition imposed in local coordinates of the element.
### Fields:
- `dofs`   -- Degrees of freedom where the boundary condition is imposed. 
- `values` -- Values imposed function. 
- `name`   -- Boundary condition label.
"""
struct LocalLoadBoundaryCondition <: AbstractLoadBoundaryCondition
    dofs::Vector{Symbol}
    values::Function
    name::Symbol
    function LocalLoadBoundaryCondition(dofs::Vector{Symbol}, values::Function, name::Symbol)
        @assert length(values(rand(1)...)) == 1 "Only a normal pressure is supported."
        new(dofs, values, name)
    end
end

"Constructor for `LocalLoadBoundaryCondition` with a string label."
LocalLoadBoundaryCondition(dofs::Vector{Symbol}, values::Function, name::String="no_labelled_bc") =
    LocalLoadBoundaryCondition(dofs, values, Symbol(name))


"Returns the dofs and the values imposed in the `GlobalLoadBoundaryCondition` `lbc` to the `AbstractFace` `f` at time `t`."
function _apply(lbc::LocalLoadBoundaryCondition, f::AbstractFace, t::Real)

    # Compute tension vector for each node 
    n = normal_direction(f)
    A = area(f)
    p = first(lbc(t)) * A * n
    num_nodes = length(nodes(f))
    p_nodal = p / num_nodes

    # Find dofs of the node corresponding to the dofs symbols of the boundary condition
    dofs_lbc = Dof[]
    [push!(dofs_lbc, dofs(n)[dof_symbol]...) for dof_symbol in dofs(lbc) for n in nodes(f)]

    # Repeat the values and build the tension vector for all dofs
    p_vec = repeat(p_nodal, outer=Int(length(dofs_lbc) / length(p_nodal)))

    @assert length(p_vec) == length(dofs_lbc)
    "The length of the tension vector must be equal to the length of the dofs vector."

    return dofs_lbc, p_vec
end

