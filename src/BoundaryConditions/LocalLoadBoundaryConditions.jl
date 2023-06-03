"""
Module to handle loads in local coordinates.
"""
module LocalLoadBoundaryConditions

using Reexport

using ..Utils
using ..BoundaryConditions
using ..Entities

@reexport import ..BoundaryConditions: apply

export Pressure

""" 
Pressure boundary condition imposed in local coordinates to the `AbstractElement` or `AbstractFace`.
This is a force along the minus normal direction of the face.
"""
Base.@kwdef struct Pressure <: AbstractNeumannBoundaryCondition
    "Degrees of freedom where the boundary condition is imposed."
    dofs::Vector{Field} = [:u]
    "Values imposed function."
    values::Function
    "Boundary condition label."
    name::Label = NO_LABEL
    # Check values is real Function
    function Pressure(dofs::Vector{Field}, values::Function, name::Label)
        @assert values(rand()) isa Real
        new(dofs, values, name)
    end
end

"Return the dofs and the values imposed in the `GlobalLoadBoundaryCondition` `lbc` to the `AbstractFace` `f` at time `t`."
function apply(lbc::Pressure, f::AbstractFace, t::Real)

    # Extract the normal vector and the area of the face
    n = normal_direction(f)
    A = area(f)

    # Compute the normal tension in global coordinates
    p = lbc(t) * (-n) * A
    num_nodes = length(nodes(f))
    p_nodal = p / num_nodes

    # Find dofs of the node corresponding to the dofs symbols of the boundary condition
    dofs_lbc = Dof[]
    [push!(dofs_lbc, dofs(n)[dof_symbol]...) for dof_symbol in dofs(lbc) for n in nodes(f)]

    # Repeat the values and build the tension vector for all dofs
    p_vec = repeat(p_nodal; outer=Int(length(dofs_lbc) / length(p_nodal)))

    @assert length(p_vec) == length(dofs_lbc)
    "The length of the tension vector must be equal to the length of the dofs vector."

    dofs_lbc, p_vec
end

end
