"""
Module defining meshes entities interface.
"""
module Mesh

using Reexport: @reexport
using StaticArrays: SVector

@reexport import ..Utils: index

export Node, coordinates, coordinates_eltype, dimension

# ======================
# Abstract Node
# ======================

abstract type AbstractNode end

""" Abstract supertype for all nodes.

An `AbstractNode` object is a point in space.

**Common methods:**

* [`coordinates`](@ref)
* [`coordinates_eltype`](@ref)
* [`index`](@ref)
* [`dimension`](@ref)
"""

const DEFAULT_NODE_INDEX = 0

"""
    Node{dim, T}
A `Node` is a point in space.
### Fields:
- `x`     -- stores the coordinates
- `index` -- stores the node index in the mesh.
"""
struct Node{dim,T} <: AbstractNode
    x::AbstractArray{T}
    index::Int
    function Node(x::AbstractArray{T}, index::Int=DEFAULT_NODE_INDEX) where {T}
        dim = length(x)
        new{dim,T}(x, index)
    end
end

Node(x::NTuple{dim,T}, index::Int=DEFAULT_NODE_INDEX) where {dim,T} = Node(SVector(x), index)

"Returns coordinates of an entity."
coordinates(n::Node) = n.x

"Returns coordinate's type of an entity."
coordinates_eltype(::Node{dim,T}) where {dim,T} = T

index(n::Node) = n.index

"Returns dimension."
dimension(::Node{dim}) where {dim} = dim

Base.getindex(n::Node, i::Int) = n.x[i]

#=

# ======================
# Abstract Mesh
# ======================

""" Abstract supertype for all meshes.

The following methods are provided by the interface:
- `element_types`    -- returns the element types that are present in the mesh. 
- `dimension`        -- returns the dimension of the mesh (1D, 2D or 3D). 
- `nodes_coordinates` -- returns a matrix with the nodes coordinates. 
- `connectivity`      -- returns the mesh connectivity. 
- `node_type (m)`    -- returns the node coordinates type.
"""

abstract type AbstractMesh{D,E,T} end

const ERROR_MESH = :("This method is not available for this mesh type. Please implement it")

" Returns the nodes data type "
node_type(::AbstractMesh{D,E,T}) = T

" Returns the mesh dimension "
dimension(::AbstractMesh{D}) where {D} = D

" Returns the element types "
element_types(::AbstractMesh{D,E}) where {D,E} = E

"Returns nodes coordinates matrix"
nodal_coordinates(::AbstractMesh) = error(ERROR_MESH)

"Returns the mesh connectivity"
connectivity(::AbstractMesh) = error(ERROR_MESH)

""" Mesh.
### Fields:
- `nodal_coords` -- Nodal coordinates matrix.
- `elements` -- Set of elements used .
"""
struct Mesh{D,E,T,C} <: AbstractMesh{D,E,T} where {D,E,T}
    nodal_coordinates::Tuple{NTuple{D,T}}
    elements::E
    connectivity::C
    MGBI_mat::Matrix{Int64}
    MGBI_vec::Vector{Int64}
end

nodal_coordinates(m::Mesh{D,E,T,C}) = m.nodal_coordinates
connectivity(m::Mesh{D,E,T,C}) = m.connectivity


=#

end # module