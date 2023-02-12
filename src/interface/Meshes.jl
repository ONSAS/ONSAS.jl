"""
Module defining meshes entities interface.
"""
module Meshes

using Reexport: @reexport
@reexport using Dictionaries: Dictionary, dictionary

using ..Utils: row_vector
@reexport using ..Elements
import ..Elements: add_dofs!
@reexport import ..Utils: dimension, dofs, nodes

export Mesh, dimension, dofs, num_dofs, elements, num_elements, element_sets, nodes, num_nodes, node_sets

""" Abstract supertype for all meshes.

The following methods are provided by the interface:


**Common methods:**

* [`dimension`](@ref)
* [`dofs`](@ref)
* [`num_dofs`](@ref)
* [`elements`](@ref)
* [`num_elements`](@ref)
* [`element_sets`](@ref)
* [`nodes`](@ref)
* [`num_nodes`](@ref)
* [`node_sets`](@ref)
"""

abstract type AbstractMesh{dim} end

"Returns the dimension of the mesh. "
dimension(::AbstractMesh{dim}) where {dim} = dim

"Returns a vector of dofs, entry `i` represents the dofs of node index `i` "
dofs(m::AbstractMesh) = dofs.(nodes(m))

"Retunrs ture if the mesh has dofs defined"
_have_dofs(m::AbstractMesh) = !all(isempty.(dofs(m)))

"Returns the maximum dof index. This function assumes that dofs start at `Dof(1)`"
function num_dofs(m::AbstractMesh)::Int
    mesh_dofs = dofs(m)
    max_dof = 0
    !_have_dofs(m) && return max_dof # mesh has no dofs

    for node_dof in mesh_dofs
        for dofs in values(node_dof)
            max_dof = maximum([max_dof, maximum(dofs)])
        end
    end
    return max_dof
end

"Add dofs to the mesh with a symbol, and considering a certain number of dofs per node."
function add_dofs!(m::AbstractMesh, dof_symbol::Symbol, dofs_per_node::Int)

    mesh_dofs = dofs(m)
    dof_not_added_yet = dof_symbol ∉ keys.(mesh_dofs)
    @assert dof_not_added_yet throw(ArgumentError("Dof symbol $dof_symbol already exists."))

    if !_have_dofs(m) && dof_not_added_yet  # any dof has been added

        for (i, n) in enumerate(nodes(m))
            node_dofs_int = (1+(i-1)*dofs_per_node):(i*dofs_per_node)
            add_dofs!(n, dof_symbol, Dof.(node_dofs_int))
        end
    else # other dof has been added
        # Maximum dof index among all dofs
        max_dof_index = num_dofs(m)
        # Push new dofs
        for (i, n) in enumerate(nodes(m))
            node_dofs_int = (1+max_dof_index+(i-1)*dofs_per_node):(max_dof_index+i*dofs_per_node)
            add_dofs!(n, dof_symbol, Dof.(node_dofs_int))
        end
    end
end

"Returns the nodes of the mesh."
nodes(m::AbstractMesh) = m.nodes

"Returns the number of nodes of the mesh."
num_nodes(m::AbstractMesh) = length(nodes(m))

"Returns the node sets of the mesh."
node_sets(m::AbstractMesh) = m.node_sets

"Returns the elements of the mesh."
elements(m::AbstractMesh) = m.elements

"Returns the number of elements of the mesh."
num_elements(m::AbstractMesh) = length(elements(m))

"Returns the element sets of the mesh."
element_sets(m::AbstractMesh) = m.element_sets

"Adds a new `Element` to the mesh"
Base.push!(m::AbstractMesh, e::AbstractElement) = push!(elements(m), e)
Base.push!(m::AbstractMesh, ve::Vector{<:AbstractElement}) = [push!(elements(m), e) for e in ve]

"Adds a new `Node` to the mesh"
Base.push!(m::AbstractMesh, n::AbstractNode) = push!(nodes(m), n)
Base.push!(m::AbstractMesh, vn::Vector{<:AbstractNode}) = [push!(nodes(m), n) for n in vn]


""" Mesh.
A `Mesh` is a collection of `Element`s and `Node`s which covers the computational domain, 
together with Sets of elements and nodes. These entities are gathered in the `element_sets` and `node_sets`.

### Fields:
- `nodes`         -- Stores the `dim` dimensional nodes of the grid.
- `elements`      -- Stores the `Element`s of the mesh.
- `node_sets      -- Maps a `String` key to a `Set` of nodes indexes.
- `element_sets   -- Maps a `String` key to a `Set` of elements indexes.

"""
struct Mesh{dim,E<:AbstractElement,N<:AbstractNode{dim}} <: AbstractMesh{dim}
    nodes::Vector{N}
    elements::Vector{E}
    # Sets
    node_sets::Dictionary{String,Set{Int}}
    element_sets::Dictionary{String,Set{Int}}
    function Mesh(
        nodes::Vector{N}=Vector{AbstractNode}(),
        elements::Vector{E}=Vector{AbstractElement}(),
        node_sets::Dictionary{String,Set{Int}}=Dictionary{String,Set{Int}}(),
        element_sets::Dictionary{String,Set{Int}}=Dictionary{String,Set{Int}}()
    ) where {N<:AbstractNode,E<:AbstractElement}

        # Check nodes dimension
        dims = dimension.(nodes)
        @assert allequal(dims) throw(ArgumentError("All nodes must have the same dimension."))

        # Check node sets indexes
        for set in keys(node_sets)
            @assert all(1 <= i <= length(nodes) for i in node_sets[set]) throw(ArgumentError("Node set $set contains invalid node indexes."))
        end

        # Check element sets indexes
        for set in keys(element_sets)
            @assert all(1 <= i <= length(elements) for i in element_sets[set]) throw(ArgumentError("Element set $set contains invalid node indexes."))
        end

        return new{first(dims),E,N}(nodes, elements, node_sets, element_sets)
    end
end



end # module