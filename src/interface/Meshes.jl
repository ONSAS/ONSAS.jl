"""
Module defining meshes entities interface.
"""
module Meshes

using Dictionaries: Dictionary
using Reexport: @reexport
using ..Utils: row_vector
using ..Elements: AbstractElement, AbstractNode, AbstractIndex, NodeIndex, ElementIndex

@reexport import ..BoundaryConditions
@reexport import ..Elements: nodes
@reexport import ..Utils: dimension, dofs

export Mesh, dimension, elements, element_nodes, element_sets, nodes, node_sets, num_dofs, num_elements

# =============
# Abstract Mesh
# =============

""" Abstract supertype for all meshes.

The following methods are provided by the interface:


**Common methods:**

* [`Base.push!`](@ref)
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

"Returns the mesh dofs "
dofs(m::AbstractMesh) = row_vector(dofs.(nodes(m)))

"Returns the number of dofs"
num_dofs(m::AbstractMesh) = length(dofs(m))

"Returns the elements of the mesh."
elements(m::AbstractMesh) = m.elements

"Returns the number of elements in the mesh"
num_elements(m::AbstractMesh) = length(elements(m))

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
Base.getindex(m::AbstractMesh, i_n::NodeIndex) = nodes(m)[i_n[]]
Base.getindex(m::AbstractMesh, vi::Vector{<:AbstractIndex}) = [m[i] for i in vi]
Base.getindex(m::AbstractMesh, si::Set{<:AbstractIndex}) = [m[i] for i in collect(si)]

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
    elements::Vector{E}
    # Sets
    node_sets::Dictionary{String,Set{NodeIndex}}
    element_sets::Dictionary{String,Set{ElementIndex}}
end

"Mesh construct with empty dofs"
function Mesh(
    nodes::Vector{N},
    elements::Vector{E}=Vector{AbstractElement}(),
    node_sets::Dictionary{String,Set{NodeIndex}}=Dictionary{String,Set{NodeIndex}}(),
    element_sets::Dictionary{String,Set{ElementIndex}}=Dictionary{String,Set{ElementIndex}}()) where {N,E}
    n₁ = first(nodes)
    D = dimension(n₁)
    T = eltype(coordinates(n₁))
    # Add nodes
    for n in enumerate(nodes)
        @assert dimension(n) != D throw(ArgumentError("All nodes must have the same dimension."))
    end
    return Mesh{D,E,N,T}(nodes, elements, node_sets, element_sets)
end

end # module