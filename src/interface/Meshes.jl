"""
Module defining meshes entities interface.
"""
module Meshes

using Reexport: @reexport
using ..Utils: NodeIndex, ElementIndex, row_vector
using ..Elements: AbstractElement, AbstractNode, set_index!


import ..Elements: coordinates_eltype
import ..Utils: dimension, dofs, elements, nodes

export Mesh, node_sets, element_nodes, element_sets, element_types, elements, connectivity

# ======================
# Abstract Mesh
# ======================

""" Abstract supertype for all meshes.

The following methods are provided by the interface:


**Common methods:**

* [`Base.push!`](@ref)
* [`coordinates_eltype`](@ref)
* [`dimension`](@ref)
* [`dofs`](@ref)
* [`elements`](@ref)
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

"Returns the mesh dofs "
dofs(m::AbstractMesh) = row_vector(dofs.(nodes(m)))

"Returns the elements of the mesh."
elements(m::AbstractMesh) = m.elements

"Returns the elements of the mesh."
element_nodes(m::AbstractMesh) = m.element_nodes

"Returns the element sets of the mesh."
element_sets(m::AbstractMesh) = m.element_sets

"Returns the nodes of the mesh."
nodes(m::AbstractMesh) = m.nodes

"Returns the node sets of the mesh."
node_sets(m::AbstractMesh) = m.node_sets

"Adds a new `Element` to the mesh"
Base.push!(m::AbstractMesh, e::AbstractElement) = push!(elements(m), e)

"Adds a new `Node` to the mesh"
Base.push!(m::AbstractMesh, n::AbstractNode) = push!(nodes(m), n)

Base.getindex(m::AbstractMesh, i_e::ElementIndex) = elements(m)[i_e[]]
Base.getindex(m::AbstractMesh, i_e::NodeIndex) = nodes(m)[i_e[]]

""" Mesh.
A `Mesh` is a collection of `Element`s and `Node`s which covers the computational domain, 
together with Sets of elements and nodes.There are multiple helper structures to apply 
boundary conditions or define subdomains. They are gathered in the `element_sets` and `node_sets`.

### Fields:
- `nodes`         -- Stores the `dim` dimensional nodes of the grid.
- `element_nodes` -- Stores the nodes of each element.
- `elements`      -- Stores the `Element`s of the mesh.
- `node_sets      -- maps a `String` key to a `Set` of nodes ids.
- `element_sets   -- maps a `String` key to a `Set` of elements ids.

"""
struct Mesh{D,E<:AbstractElement,N<:AbstractNode,T} <: AbstractMesh{D,T}
    nodes::Vector{N}
    element_nodes::Vector{Vector{N}}
    elements::Vector{E}
    # Sets
    node_sets::Dict{String,Set{NodeIndex}}
    element_sets::Dict{String,Set{ElementIndex}}
    function Mesh(
        vnodes::Vector{N},
        element_nodes::Vector{Vector{N}},
        elements::Vector{E}=Vector{AbstractElement}(),
        node_sets::Dict{String,Set{NodeIndex}}=Dict{String,Set{NodeIndex}}(),
        element_sets::Dict{String,Set{ElementIndex}}=Dict{String,Set{ElementIndex}}()) where {N,E}
        n₁ = first(vnodes)
        D = dimension(n₁)
        T = coordinates_eltype(n₁)
        # Add nodes
        for (i, n) in enumerate(vnodes)
            set_index!(n, i)
            dimension(n) != D && throw(ArgumentError("All nodes must have the same dimension."))
            coordinates_eltype(n) != T && throw(ArgumentError("All nodes must have the same eltype."))
        end
        # Add Elements
        for (i, e) in enumerate(elements)
            set_index!(e, i)
            # Check nodes ∈ nodes
            [n ∈ vnodes || throw(ArgumentError("The element node ∉ nodes(mesh)")) for n in nodes(elements)]
        end
        return new{D,E,N,T}(vnodes, element_nodes, elements, node_sets, element_sets)
    end
end

end # module