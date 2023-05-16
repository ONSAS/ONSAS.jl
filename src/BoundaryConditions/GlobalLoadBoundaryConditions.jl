"""
Module to handle loads in global coordinates.
"""
module GlobalLoadBoundaryConditions

using Reexport

using ..Utils
using ..BoundaryConditions
using ..Elements

@reexport import ..BoundaryConditions: apply

export GlobalLoad

""" 
Load boundary condition imposed in global coordinates of the element.
"""
Base.@kwdef struct GlobalLoad <: AbstractNeumannBoundaryCondition
    "Degrees of freedom where the boundary condition is imposed."
    dofs::Vector{Field} = [:u]
    "Values imposed function."
    values::Function
    "Label of the boundary condition."
    name::Label = NO_LABEL
end

"Return the dofs and the values imposed in the `GlobalLoad` `lbc` to 
the `AbstractNode` `n` at time `t`. "
function apply(lbc::GlobalLoad, n::AbstractNode, t::Real)

    # Find dofs of the node corresponding to the dofs symbols of the boundary condition
    dofs_lbc = Dof[]
    [push!(dofs_lbc, dofs(n)[dof_symbol]...) for dof_symbol in dofs(lbc)]
    f = lbc(t)

    # Repeat the values and build the force vector for all dofs
    f_dofs = repeat(f; outer=Int(length(dofs_lbc) / length(f)))

    @assert length(f_dofs) == length(dofs_lbc)
    "The length of the force vector must be equal to the length of the dofs vector."

    dofs_lbc, f_dofs
end

"Return the dofs and the values imposed in the `GlobalLoad` `lbc` to the `AbstractFace` `f` at time `t`."
function apply(lbc::GlobalLoad, f::AbstractFace, t::Real)

    # Compute tension vector for each node 
    A = area(f)
    p = lbc(t) * A
    num_nodes = length(nodes(f))
    p_nodal = p / num_nodes

    # Find dofs of the node corresponding to the dofs symbols of the boundary condition
    dofs_lbc = Dof[]
    [push!(dofs_lbc, dofs(n)[dof_symbol]...) for dof_symbol in dofs(lbc) for n in nodes(f)]

    # Repeat the values and build the tension vector for all dofs
    p_vec = repeat(p_nodal; outer=Int(length(dofs_lbc) / length(p)))

    @assert length(p_vec) == length(dofs_lbc)
    "The length of the tension vector must be equal to the length of the dofs vector."

    dofs_lbc, p_vec
end

"Return the dofs and the values imposed in the `GlobalLoad` `lbc` to the `AbstractElement` `e` at time `t`."
function apply(lbc::GlobalLoad, e::AbstractElement, t::Real)

    # Compute tension vector for each node 
    V = volume(e)
    b = lbc(t) * V
    num_nodes = length(nodes(e))
    b_nodal = b / num_nodes

    # Find dofs of the node corresponding to the dofs symbols of the boundary condition
    dofs_lbc = Dof[]
    [push!(dofs_lbc, dofs(n)[dof_symbol]...) for dof_symbol in dofs(lbc) for n in nodes(e)]

    # Repeat the values and build the tension vector for all dofs
    b_vec = repeat(b_nodal; outer=Int(length(dofs_lbc) / length(b_nodal)))

    @assert length(b_vec) == length(dofs_lbc)
    "The length of the tension vector must be equal to the length of the dofs vector."

    dofs_lbc, b_vec
end

end
