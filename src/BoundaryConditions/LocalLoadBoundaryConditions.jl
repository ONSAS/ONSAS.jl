"""
Module to handle loads in local coordinates.
"""
module LocalLoadBoundaryConditions

using Reexport

using ..Utils
using ..BoundaryConditions
using ..Entities
using ..Nodes

@reexport import ..BoundaryConditions: apply

export Pressure

"""
Pressure boundary condition imposed in local coordinates to the `AbstractElement` or `AbstractFace`.
This is a force along the minus normal direction of the face.
"""
struct Pressure <: AbstractNeumannBoundaryCondition
    "Field where the boundary condition applies to."
    field::Field
    "Values imposed function."
    values::Function
    "Boundary condition label."
    name::Label
    function Pressure(field::Field, values::Function, name::Label=NO_LABEL)
        # Check values is real function.
        @assert values(rand()) isa Real
        new(field, values, name)
    end
end

"Return the dofs and the values imposed in the `GlobalLoadBoundaryCondition` `lbc` to the `AbstractFace` `f` at time `t`."
function apply(bc::Pressure, f::AbstractFace, t::Real)

    # Extract the normal vector and the area of the face
    n = normal_direction(f)
    A = area(f)

    # Compute the normal tension in global coordinates
    p = bc(t) * (-n) * A
    num_nodes = length(nodes(f))
    p_nodal = p / num_nodes

    # Find dofs of the node corresponding to the dofs symbols of the boundary condition
    dofs_bc = reduce(vcat, values(dofs(n)[bc.field]) for n in nodes(f))

    # Repeat the values and build the tension vector for all dofs
    p_vec = repeat(p_nodal; outer=Int(length(dofs_bc) / length(p_nodal)))

    @assert length(p_vec) == length(dofs_bc) "The length of the tension vector must be equal to the length of the dofs vector."

    dofs_bc, p_vec
end

end
