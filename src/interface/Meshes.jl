"""
Module defining meshes entities interface.
"""
module Meshes

using Reexport: @reexport

export Mesh, build_node_set, nodes, node_sets, elements,
    element_sets, element_types, elements, connectivity

# ======================
# Abstract Mesh
# ======================

""" Abstract supertype for all meshes.

The following methods are provided by the interface:


**Common methods:**

* [`connectivity_matrix`](@ref)
* [`coordinates`](@ref)
* [`coordinates_eltype`](@ref)
* [`elements`](@ref)
* [`element_types`](@ref)
* [`element_sets`](@ref)
* [`nodes`](@ref)
* [`node_sets`](@ref)
* [`dimension`](@ref)
"""

abstract type AbstractMesh{D,E,T} end

"Returns the dimension of the mesh. "
dimension(::AbstractMesh{D}) where {D} = D

"Returns the element types in the mesh."
element_types(::AbstractMesh{D,E}) where {D,E} = E

"Returns the coordinate's type of the mesh. "
coordinates_eltype(::AbstractMesh{D,E,T}) where {D,E,T} = T


""" Mesh.
A `Mesh` is a collection of `Elements` and `Node`s which covers the computational domain, together with Sets of elements and nodes.
There are multiple helper structures to apply boundary conditions or define subdomains. They are gathered in the `elementsets` and `nodesets`.

### Fields:
- `nodes`       -- Stores the `dim` dimensional nodes of the grid
- `elements`    -- Stores all elements of the grid
- `node_sets    -- maps a `String` key to a `Set` of nodes ids.
- `element_sets -- maps a `String` key to a `Set` of elements ids.

"""
Base.@kwdef struct Mesh{D,E,N,T} <: AbstractMesh{D,E,T}
    nodes::Vector{N}
    elements::Vector{E}
    # Sets
    node_sets::Dict{String,Set{Int}} = Dict{String,Set{Int}}()
    element_sets::Dict{String,Set{Int}} = Dict{String,Set{Int}}()
    function Mesh(
        nodes::Vector{N},
        elements::Vector{E},
        node_sets::Dict{String,Set{Int}}=Dict{String,Set{Int}}(),
        element_sets::Dict{String,Set{Int}}=Dict{String,Set{Int}}()) where {E,N}
        n₁ = first(nodes)
        D = dimension(n₁)
        T = coordinates_eltype(n₁)
        # Assign indexes
        for (i, n) in enumerate(nodes)
            set_index!(n, i)
            dimension(n) != D && throw(ArgumentError("All nodes must have the same dimension."))
            coordinates_eltype(n) != T && throw(ArgumentError("All nodes must have the same eltype."))
        end
        # [set_index!(e, i) for (i, e) in enumerate(elements)]
        return new{D,E,N,T}(nodes, elements, node_sets, element_sets)
    end
end

"Returns the nodes of the mesh."
nodes(m::Mesh) = m.nodes

"Returns the node sets of the mesh."
node_sets(m::Mesh) = m.node_sets

"Create a nodes set in a Tuple (label, set)"
build_node_set(label::String, n::Vector{N}) where {N} = (label, Set(index.(n)))

"Returns the elements of the mesh."
elements(m::Mesh) = m.elements

"Returns the element sets of the mesh."
element_sets(m::Mesh) = m.element_sets

"Create a element set in a Tuple (label, set)"
build_elements_set(label::String, e::Vector{E}) where {E} = (label, Set(index.(e)))


"Connectivity matrix"
function connectivity_matrix(m::Mesh{D,E,T}) where {D,E,T}

    elems = elements(m)
    C = zeros(Int, length(elems), length(m.elements[1].nodes))
    for e in elements(m)

        # TO BE DONE
    end
    return C
end


end # module