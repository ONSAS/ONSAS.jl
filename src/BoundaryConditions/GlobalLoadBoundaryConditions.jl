"""
Module to handle loads in global coordinates.
"""
module GlobalLoadBoundaryConditions

using Reexport

using ..Utils
using ..BoundaryConditions
using ..Nodes
using ..TriangularFaces
using ..Entities

@reexport import ..BoundaryConditions: apply

export GlobalLoad

"""
Load boundary condition imposed in global coordinates of the element.
"""
struct GlobalLoad <: AbstractNeumannBoundaryCondition
    "Field where the boundary condition applies to."
    field::Field
    "Values imposed function."
    values::Function
    "Boundary condition label."
    name::Label
    function GlobalLoad(field::Field, values::Function, name::Label=NO_LABEL)
        new(field, values, name)
    end
end

"Return the dofs and the values imposed in the global load to the `AbstractNode` `n` at time `t`. "
function apply(bc::GlobalLoad, n::AbstractNode, t::Real)

    # Find dofs of the node corresponding to the dofs symbols of the boundary condition
    dofs_bc = dofs(n)[bc.field]
    f = bc(t)

    # Repeat the values and build the force vector for all dofs
    f_dofs = repeat(f; outer=Int(length(dofs_bc) / length(f)))

    @assert length(f_dofs) == length(dofs_bc)
    "The length of the force vector must be equal to the length of the dofs vector."

    dofs_bc, f_dofs
end

"Return the dofs and the values imposed in the global load to the `AbstractFace` `f` at time `t`."
function apply(bc::GlobalLoad, f::AbstractFace, t::Real)

    # Compute tension vector for each node
    A = area(f)
    p = bc(t) * A
    num_nodes = length(nodes(f))
    p_nodal = p / num_nodes

    # Find dofs of the node corresponding to the dofs symbols of the boundary condition
    dofs_bc = reduce(vcat, values(dofs(n)[bc.field]) for n in nodes(f))

    # Repeat the values and build the tension vector for all dofs
    p_vec = repeat(p_nodal; outer=Int(length(dofs_bc) / length(p)))

    @assert length(p_vec) == length(dofs_bc)
    "The length of the tension vector must be equal to the length of the dofs vector."

    dofs_bc, p_vec
end

"Return the dofs and the values imposed in the global load to the `AbstractElement` `e` at time `t`."
function apply(bc::GlobalLoad, e::AbstractElement, t::Real)

    # Compute tension vector for each node
    V = volume(e)
    b = bc(t) * V
    num_nodes = length(nodes(e))
    b_nodal = b / num_nodes

    # Find dofs of the node corresponding to the dofs symbols of the boundary condition
    dofs_bc = reduce(vcat, values(dofs(n)[bc.field]) for n in nodes(e))

    # Repeat the values and build the tension vector for all dofs
    b_vec = repeat(b_nodal; outer=Int(length(dofs_bc) / length(b_nodal)))

    @assert length(b_vec) == length(dofs_bc)
    "The length of the tension vector must be equal to the length of the dofs vector."

    dofs_bc, b_vec
end

end
