using ..BoundaryConditions: AbstractLoadBoundaryCondition
using ..Elements: AbstractNode, AbstractFace, AbstractElement, area, volume

import ..BoundaryConditions: _apply

export GlobalLoadBoundaryCondition

""" Load boundary condition imposed in global coordinates of the element.
### Fields:
- `dofs`   -- Degrees of freedom where the boundary condition is imposed. 
- `values` -- Values imposed function. 
- `name`   -- Boundary condition label.
"""
struct GlobalLoadBoundaryCondition <: AbstractLoadBoundaryCondition
    dofs::Vector{Symbol}
    values::Function
    name::Symbol
end

"Constructor for `LocalLoadBoundaryCondition` with a string label."
GlobalLoadBoundaryCondition(dofs::Vector{Symbol}, values::Function, name::String="no_labelled_bc") =
    GlobalLoadBoundaryCondition(dofs, values, Symbol(name))

"Returns the dofs and the values imposed in the `GlobalLoadBoundaryCondition` `lbc` to 
the `AbstractNode` `n` at time `t`. "
function _apply(lbc::GlobalLoadBoundaryCondition, n::AbstractNode, t::Real)

    # Find dofs of the node corresponding to the dofs symbols of the boundary condition
    dofs_lbc = Dof[]
    [push!(dofs_lbc, dofs(n)[dof_symbol]...) for dof_symbol in dofs(lbc)]
    f = lbc(t)

    # Repeat the values and build the force vector for all dofs
    f_dofs = repeat(f, outer=Int(length(dofs_lbc) / length(f)))

    @assert length(f_dofs) == length(dofs_lbc)
    "The length of the force vector must be equal to the length of the dofs vector."

    return dofs_lbc, f_dofs
end

"Returns the dofs and the values imposed in the `GlobalLoadBoundaryCondition` `lbc` to the `AbstractFace` `f` at time `t`."
function _apply(lbc::GlobalLoadBoundaryCondition, f::AbstractFace, t::Real)

    # Compute tension vector for each node 
    A = area(f)
    p = lbc(t) * A
    num_nodes = length(nodes(f))
    p_nodal = p / num_nodes

    # Find dofs of the node corresponding to the dofs symbols of the boundary condition
    dofs_lbc = Dof[]
    [push!(dofs_lbc, dofs(n)[dof_symbol]...) for dof_symbol in dofs(lbc) for n in nodes(f)]

    # Repeat the values and build the tension vector for all dofs
    p_vec = repeat(p_nodal, outer=Int(length(dofs_lbc) / length(p)))

    @assert length(p_vec) == length(dofs_lbc)
    "The length of the tension vector must be equal to the length of the dofs vector."

    return dofs_lbc, p_vec
end


"Returns the dofs and the values imposed in the `GlobalLoadBoundaryCondition` `lbc` to the `AbstractElement` `e` at time `t`."
function _apply(lbc::GlobalLoadBoundaryCondition, e::AbstractElement, t::Real)

    # Compute tension vector for each node 
    V = volume(e)
    b = lbc(t) * V
    num_nodes = length(nodes(e))
    b_nodal = b / num_nodes

    # Find dofs of the node corresponding to the dofs symbols of the boundary condition
    dofs_lbc = Dof[]
    [push!(dofs_lbc, dofs(n)[dof_symbol]...) for dof_symbol in dofs(lbc) for n in nodes(e)]

    # Repeat the values and build the tension vector for all dofs
    b_vec = repeat(b_nodal, outer=Int(length(dofs_lbc) / length(b_nodal)))

    @assert length(b_vec) == length(dofs_lbc)
    "The length of the tension vector must be equal to the length of the dofs vector."

    return dofs_lbc, b_vec
end

