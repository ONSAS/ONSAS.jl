"""
Module defining meshes entities interface.
Each mesh consists of a data type with the nodes and elements. Moreover, different sets of nodes and elements can be defined.
"""
module Meshes

using Reexport: @reexport

@reexport using ..Elements
using ..Utils: row_vector

@reexport import ..Elements: add!, dimension, dofs, nodes

export AbstractMesh, Mesh, elements, element_sets, num_dofs, num_elements, num_nodes, node_sets

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

"Returns the dimension of an `AbstractMesh`."
dimension(::AbstractMesh{dim}) where {dim} = dim

"Returns the `AbstractMesh` vector of `Dof`s. 
Entry `i` contains the `Dof`s of node with index `i` in the `AbstractMesh` vector of nodes."
dofs(m::AbstractMesh) = dofs.(nodes(m))

"Returns true if the `AbstractMesh` m has `Dof`s defined."
_isempty_dofs(m::AbstractMesh) = !all(isempty.(dofs(m)))

"Returns the number of `Dof`s defined in the `AbstractMesh` `m`.
This function assumes that `Dof`s indexes start from `Dof(1)`"
function num_dofs(m::AbstractMesh)::Int
    mesh_dofs = dofs(m)
    max_dof = 0
    !_isempty_dofs(m) && return max_dof # mesh has no dofs

    for node_dof in mesh_dofs
        for dofs in values(node_dof)
            max_dof = maximum([max_dof, maximum(dofs)])
        end
    end
    return max_dof
end

"Adds n `dofs_per_node` `Dof`s with `dof_symbol` to the `AbstractMesh` `m`."
function add!(m::AbstractMesh, dof_symbol::Symbol, dofs_per_node::Int)

    mesh_dofs = dofs(m)
    dof_not_added_yet = dof_symbol ∉ keys.(mesh_dofs)
    @assert dof_not_added_yet throw(ArgumentError("Dof symbol $dof_symbol already exists."))

    if !_isempty_dofs(m) && dof_not_added_yet  # any dof has been added

        for (i, n) in enumerate(nodes(m))
            node_dofs_int = (1+(i-1)*dofs_per_node):(i*dofs_per_node)
            add!(n, dof_symbol, Dof.(node_dofs_int))
        end
    else # other dof has been added
        # Maximum dof index among all dofs
        max_dof_index = num_dofs(m)
        # Push new dofs
        for (i, n) in enumerate(nodes(m))
            node_dofs_int = (1+max_dof_index+(i-1)*dofs_per_node):(max_dof_index+i*dofs_per_node)
            add!(n, dof_symbol, Dof.(node_dofs_int))
        end
    end
end

"Returns the `Node`s of the `AbstractMesh` `m`."
nodes(m::AbstractMesh) = m.nodes

"Returns the number of `Node`s of the `AbstractMesh` `m`."
num_nodes(m::AbstractMesh) = length(nodes(m))

"Returns `Node` sets of the `AbstractMesh` `m`."
node_sets(m::AbstractMesh) = m.node_sets

"Returns the `Element`s of the `AbstractMesh` `m`."
elements(m::AbstractMesh) = m.elements

"Returns the number of elements of the `AbstractMesh` `m`."
num_elements(m::AbstractMesh) = length(elements(m))

"Returns the element sets of the `AbstractMesh` `m`."
element_sets(m::AbstractMesh) = m.element_sets

"Pushes a new `Element` into the `AbstractMesh` `m`."
Base.push!(m::AbstractMesh, e::AbstractElement) = push!(elements(m), e)

"Pushes a new vector of `Element`s into the `AbstractMesh` `m`."
Base.push!(m::AbstractMesh, ve::Vector{<:AbstractElement}) = [push!(elements(m), e) for e in ve]

"Pushes a new `Node` into the `AbstractMesh` `m`."
Base.push!(m::AbstractMesh, n::AbstractNode) = push!(nodes(m), n)

"Pushes a new  vector of `Node`s into the `AbstractMesh` `m`."
Base.push!(m::AbstractMesh, vn::Vector{<:AbstractNode}) = [push!(nodes(m), n) for n in vn]


""" Mesh.
A `Mesh` is a collection of `Element`s and `Node`s that cover the discretized domain, 
together with Sets of elements and nodes. These entities are gathered in the `element_sets` and `node_sets`.

### Fields:
- `nodes`         -- Stores the `dim` dimensional nodes of the grid.
- `elements`      -- Stores the `Element`s of the mesh.
- `node_sets`     -- Maps a `String` key to a `Set` of nodes indexes.
- `element_sets`   -- Maps a `String` key to a `Set` of elements indexes.
"""
struct Mesh{dim,E<:AbstractElement,N<:AbstractNode{dim}} <: AbstractMesh{dim}
    # Entities
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
            check = all([i ≤ length(nodes) for i in node_sets[set]])
            if check == false
                throw(ArgumentError("Node set: $set contains invalid node indexes."))
            end
        end

        # Check element sets indexes
        for set in keys(element_sets)
            check = all([i ≤ length(elements) for i in element_sets[set]])
            if check == false
                throw(ArgumentError("Element set: $set contains invalid element indexes."))
            end
        end

        return new{first(dims),E,N}(nodes, elements, node_sets, element_sets)

    end
end


end # module