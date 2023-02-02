"""
Module defining meshes entities interface.
"""
module Meshes

using Reexport: @reexport
using ..Elements: set_index!


import ..Elements: coordinates_eltype
import ..Utils: dimension, nodes

export Mesh, build_node_set, node_sets, element_nodes,
    element_sets, element_types, elements, connectivity

# ======================
# Abstract Mesh
# ======================

""" Abstract supertype for all meshes.

The following methods are provided by the interface:


**Common methods:**

* [`coordinates_eltype`](@ref)
* [`dimension`](@ref)
* [`element_nodes`](@ref)
* [`element_sets`](@ref)
* [`nodes`](@ref)
* [`node_sets`](@ref)
"""

abstract type AbstractMesh{D,T} end

"Returns the dimension of the mesh. "
dimension(::AbstractMesh{D}) where {D} = D

"Returns the coordinate's type of the mesh. "
coordinates_eltype(::AbstractMesh{D,T}) where {D,T} = T


""" Mesh.
A `Mesh` is a collection of `Element`s and `Node`s which covers the computational domain, 
together with Sets of elements and nodes.There are multiple helper structures to apply 
boundary conditions or define subdomains. They are gathered in the `element_sets` and `node_sets`.

### Fields:
- `nodes`         -- Stores the `dim` dimensional nodes of the grid.
- `element_nodes` -- Stores the nodes of each element.
- `node_sets      -- maps a `String` key to a `Set` of nodes ids.
- `element_sets   -- maps a `String` key to a `Set` of elements ids.

"""
struct Mesh{D,N,T} <: AbstractMesh{D,T}
    nodes::Vector{N}
    element_nodes::Vector{Vector{N}}
    # Sets
    node_sets::Dict{String,Set{Int}}
    element_sets::Dict{String,Set{Int}}
    function Mesh(
        nodes::Vector{N},
        element_nodes::Vector{Vector{N}},
        node_sets::Dict{String,Set{Int}}=Dict{String,Set{Int}}(),
        element_sets::Dict{String,Set{Int}}=Dict{String,Set{Int}}()) where {N}
        n₁ = first(nodes)
        D = dimension(n₁)
        T = coordinates_eltype(n₁)
        # Assign indexes
        for (i, n) in enumerate(nodes)
            set_index!(n, i)
            dimension(n) != D && throw(ArgumentError("All nodes must have the same dimension."))
            coordinates_eltype(n) != T && throw(ArgumentError("All nodes must have the same eltype."))
        end
        return new{D,N,T}(nodes, element_nodes, node_sets, element_sets)
    end
end

"Returns the nodes of the mesh."
nodes(m::Mesh) = m.nodes

"Returns the elements of the mesh."
element_nodes(m::Mesh) = m.element_nodes

"Returns the node sets of the mesh."
node_sets(m::Mesh) = m.node_sets

"Returns the element sets of the mesh."
element_sets(m::Mesh) = m.element_sets

end # module