"""
Module defining the elements implemented.
"""
module Elements

using Reexport: @reexport
using StaticArrays: SVector
using ..Utils: row_vector

@reexport using ..CrossSections
@reexport import ..Utils: dimension, label
@reexport import ..Utils: internal_forces, inertial_forces
@reexport import ..Utils: coordinates, dofs, index, nodes

export Dof
export AbstractNode, Node, coordinates
export AbstractElement, cross_section, coordinates, local_dofs

# ========================
# Degree of freedom (Dof)
# ========================
const _DEFAULT_INDEX_INT = 0

""" Degree of freedom struct.
This is a scalar degree of freedom of the structure.
### Fields:
- `index`   -- degree of freedom identification number. 
"""
mutable struct Dof
    index::Int
end

index(d::Dof) = d.index
Base.setindex!(d::Dof, i::Int) = d.index = i

@inline Base.getindex(v::AbstractVector, d::Dof) = v[index(d)]
@inline Base.getindex(v::AbstractVector, vd::Vector{<:Dof}) = [v[index(d)] for d in vd]
@inline Base.setindex!(v::AbstractVector{T}, t::T, d::Dof) where {T} = v[index(d)] = t
@inline Base.setindex!(v::AbstractVector, tv::Vector{T}, vd::Vector{<:Dof}) where {T} = [setindex!(v, ti, vi) for (ti, vi) in zip(tv, vd)]

# =================
# Abstract Node
# =================

abstract type AbstractNode{dim,T} end

""" Abstract supertype for all nodes.

An `AbstractNode` object is a point in space.

**Common methods:**

* [`coordinates`](@ref)
* [`dimension`](@ref)
* [`dofs`](@ref)
"""

coordinates(n::AbstractNode) = n.x

dimension(::AbstractNode{dim}) where {dim} = dim

Base.getindex(n::AbstractNode, i::Int) = n.x[i]

"Returns node's degrees of freedom."
dofs(n::AbstractNode) = n.dofs
dofs(n::Vector{<:AbstractNode}) = vcat(dofs.(n)...)

"Returns corresponding degrees of freedom considering angles"
function _nodes2dofs(i::Integer, dim::Integer)
    range = if dim == 3
        (i-1)*2dim+1:(i)*2dim
    elseif dim == 2
        (i-1)*2dim:(i)*2dim-1
    end
end

function Base.setindex!(n::AbstractNode{dim}, i::Int) where {dim}
    ndofs = dofs(n)
    node_dof_indexes = _nodes2dofs(i, dim)
    [setindex!(dof, node_dof_indexes[i]) for (i, dof) in enumerate(ndofs)]
    return n
end

#TODO: generalize to any field type and dimension (:u, dim)
"Maps dimension of the node to local degrees of freedom."
function _dim_to_nodal_dofs(dim::Int)
    dofs = if dim == 1
        [Dof(1)] #uₓ
    elseif dim == 2
        [Dof(1), Dof(2), Dof(3)]#uᵢ, uⱼ, θₖ
    elseif dim == 3
        [
            Dof(1), Dof(2),
            Dof(3), Dof(4),
            Dof(5), Dof(6)
        ]#uᵢ, θᵢ, uⱼ, θⱼ, uₖ, θₖ
    else
        error("Dimension not supported.")
    end
end

"""
A `Node` is a point in space.
### Fields:
- `x`     -- stores the coordinates.
- `dofs`  -- stores the node degrees of freedom.
"""
struct Node{dim,T} <: AbstractNode{dim,T}
    x::AbstractArray{T}
    dofs::Vector{<:Dof}
    function Node(
        x::AbstractArray{T}, dofs::Vector{<:Dof}) where {T<:Real}
        dim = length(x)
        @assert dim ≤ 3 "Only 1D,2D or 3D nodes are supported"
        new{dim,T}(x, dofs)
    end
end

Node(x₁::T) where {T<:Real} = Node(SVector(x₁), _dim_to_nodal_dofs(2))
Node(x₁::T, x₂::T) where {T<:Real} = Node(SVector((x₁, x₂)), _dim_to_nodal_dofs(2))
Node(x₁::T, x₂::T, x₃::T) where {T<:Real} = Node(SVector((x₁, x₂, x₃)), _dim_to_nodal_dofs(3))
Node(x::NTuple{dim,T}) where {dim,T<:Real} = Node(SVector(x), _dim_to_nodal_dofs(length(x)))
Node(x::Vector{T}) where {T<:Real} = Node(SVector(x...), _dim_to_nodal_dofs(length(x)))
# =================
# Abstract Element
# =================

const _DEFAULT_LABEL = :no_labelled_element

#TODO: Add interpolation order
abstract type AbstractElement{dim,T} end

""" Abstract supertype for all elements.

An `AbstractElement` object facilitates the process of evaluating:

    - The internal forces vector and its tangent matrices.
    - The inertial forces vector and its tangent matrices.
    - The mechanical stresses and strains.

**Common methods:**

* [`coordinates`](@ref)
* [`dofs`](@ref)
* [`local_dofs`](@ref)
* [`label`](@ref)
* [`nodes`](@ref)

* [`internal_forces`](@ref)
* [`inertial_forces`](@ref)

**Common fields:**
* nodes
* label

**Hard contracts:**

* [`local_dofs`](@ref)      - defines the local dofs of the element.

For static cases the following methods are required:

* [`inertial_forces`](@ref) - function that returns the internal forces vector and its respective tangent matrices.

For dynamic cases the following methods are required:
* [`inertial_forces`](@ref) - function that returns the inertial forces vector and its respective tangent matrices.
"""

"Returns element coordinates."
coordinates(e::AbstractElement) = row_vector(coordinates.(nodes(e)))

"Returns the geometrical properties of the element"
cross_section(e::AbstractElement) = e.cross_section

dofs(e::AbstractElement) = row_vector(dofs.(nodes(e)))
dofs(ve::Vector{<:AbstractElement}) = unique(row_vector(dofs.(ve)))

"Returns local dofs of an element. This dofs are essential for the assemble process."
function local_dofs(e::AbstractElement) end

index(e::AbstractElement) = e.index

"Returns element label."
label(e::AbstractElement) = e.label

"Returns element nodes."
nodes(e::AbstractElement) = e.nodes

"Returns the internal force vector of the element."
function internal_forces(e::AbstractElement, args...; kwargs...) end

"Returns the inertial force vector of the element."
function inertial_forces(e::AbstractElement, args...; kwargs...) end

"Returns the element stresses"
function stress(e::AbstractElement, args...; kwargs...) end

"Returns the element strain"
function strain(e::AbstractElement, args...; kwargs...) end

include("../elements/Truss.jl")

end # module


