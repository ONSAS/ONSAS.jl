"""
Module defining the elements implemented.
"""
module Elements

using AutoHashEquals: @auto_hash_equals
using Reexport: @reexport
using StaticArrays: SVector, SMatrix

using ..Materials: AbstractMaterial
using ..BoundaryConditions: AbstractBoundaryCondition
using ..InitialConditions: AbstractInitialCondition
using ..Utils: ScalarWrapper
@reexport import ..Utils: dimension, label, set_label!
@reexport import ..Utils: internal_forces, inertial_forces
@reexport import ..BoundaryConditions: dofs

export AbstractIndex, ElementIndex, NodeIndex, index
export Dof
export AbstractNode, Node, boundary_conditions, coordinates, coordinates_eltype
export AbstractElement, create_element, coordinates, dofs, local_dofs, material

const _DEFAULT_LABEL = :no_labelled_element

# ========
# Indexes 
# ========

""" Abstract supertype for all indexes.

The following methods are provided by the interface:


**Common methods:**

* [`Base.getindex`](@ref)
* [`Base.setindex!`](@ref)
* [`Base.isequal`](@ref)

"""

abstract type AbstractIndex{I} end

@inline Base.getindex(i::AbstractIndex) = i.id
@inline Base.setindex!(i::AbstractIndex, id) = i.id = id # callable with i[] = id

"""
`Element` identification number.
### Fields:
`id` -- integer number. 
"""
@auto_hash_equals mutable struct ElementIndex{I<:Integer} <: AbstractIndex{I}
    id::I
end

"""
`Node` identification number.
### Fields:
`id` -- index number. 
"""
@auto_hash_equals mutable struct NodeIndex{I<:Integer} <: AbstractIndex{I}
    id::I
end

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

const _DEFAULT_NODE_INDEX = NodeIndex(_DEFAULT_INDEX_INT)
const _DEFAULT_ELEMENT_INDEX = ElementIndex(_DEFAULT_INDEX_INT)
abstract type AbstractNode{dim,T} end

""" Abstract supertype for all nodes.

An `AbstractNode` object is a point in space.

**Common methods:**

* [`boundary_conditions`](@ref)
* [`coordinates`](@ref)
* [`coordinates_eltype`](@ref)
* [`dimension`](@ref)
* [`index`](@ref)
* [`dofs`](@ref)
"""

"Returns the node boundary conditions."
boundary_conditions(n::AbstractNode) = n.bc

Base.push!(n::AbstractNode, bc::AbstractBoundaryCondition) = push!(n.bc, bc)
Base.push!(n::AbstractNode, vbc::Vector{<:AbstractBoundaryCondition}) = [push!(n, bc) for bc in vbc]
Base.push!(vn::Vector{<:AbstractNode}, bc::AbstractBoundaryCondition) = [push!(n.bc, bc) for n in vn]

coordinates(n::AbstractNode) = n.x

"Returns coordinate's type of an entity."
coordinates_eltype(::AbstractNode{dim,T}) where {dim,T} = T

dimension(::AbstractNode{dim}) where {dim} = dim

Base.getindex(n::AbstractNode, i::Int) = n.x[i]

"Returns node's degrees of freedom."
dofs(n::AbstractNode) = n.dofs

dofs(n::Vector{<:AbstractNode}) = vcat(dofs.(n)...)

index(n::AbstractNode) = n.index[]

_nodes2dofs(idx::Integer, dim::Integer) = (idx-1)*2dim+1:(idx)*2dim

function Base.setindex!(n::AbstractNode{dim}, idx::NodeIndex) where {dim}
    ndofs = dofs(n)
    node_dof_indexes = _nodes2dofs(idx[], dim)
    [setindex!(dof, node_dof_indexes[i]) for (i, dof) in enumerate(ndofs)]
    n.index[] = idx[]
    return n
end

Base.setindex!(n::AbstractNode, i::Integer) = setindex!(n, NodeIndex(i))

#TODO: generalize to any field type and dimension (:u, dim)
"Maps dimension of the node to local degrees of freedom."
function _dim_to_nodal_dofs(dim::Int)
    dofs = if dim == 1
        [Dof(1)]
    elseif dim == 2
        [Dof(1), Dof(2), Dof(3)]
    elseif dim == 3
        [
            Dof(1), Dof(2),
            Dof(3), Dof(4),
            Dof(5), Dof(6)
        ]
    else
        error("Dimension not supported.")
    end
end

"""
A `Node` is a point in space.
### Fields:
- `x`     -- stores the coordinates.
- `dofs`  -- stores the node degrees of freedom.
- `index` -- stores the node index in the mesh.
- `bc`    -- stores the boundary conditions applied on the node.
- `ic`    -- stores the initial conditions applied on the node.
"""
mutable struct Node{dim,T} <: AbstractNode{dim,T}
    x::AbstractArray{T}
    dofs::Vector{<:Dof}
    index::NodeIndex
    bc::Vector{<:AbstractBoundaryCondition}
    ic::Vector{<:AbstractInitialCondition}
    function Node(
        x::AbstractArray{T}, dofs::Vector{<:Dof}, index::NodeIndex,
        bc::Vector{<:AbstractBoundaryCondition},
        ic::Vector{<:AbstractInitialCondition}) where {T}
        new{length(x),T}(x, dofs, index, bc, ic)
    end
end

Node(x::NTuple{dim,T}, index::Integer, bc::Vector{<:AbstractBoundaryCondition}) where {dim,T} = Node(SVector(x), _dim_to_nodal_dofs(length(x)), Index(index), bc)
# Empty bc and ic constructors
empty_bc = Vector{AbstractBoundaryCondition}()
empty_ic = Vector{AbstractInitialCondition}()
Node(x::NTuple{dim,T}, index::NodeIndex) where {dim,T} = Node(SVector(x), _dim_to_nodal_dofs(length(x)), index, empty_bc, empty_ic)
Node(x::NTuple{dim,T}, index::Integer=_DEFAULT_INDEX_INT) where {dim,T} = Node(SVector(x), _dim_to_nodal_dofs(length(x)), NodeIndex(index), empty_bc, empty_ic)
Node(x::AbstractArray{T}, index::NodeIndex) where {T} = Node(x, _dim_to_nodal_dofs(length(x)), index, empty_bc, empty_ic)
Node(x::AbstractArray{T}, index::Integer=_DEFAULT_INDEX_INT) where {T} = Node(x, _dim_to_nodal_dofs(length(x)), NodeIndex(index), empty_bc, empty_ic)

# =================
# Abstract Element
# =================

#TODO: Add interpolation order
abstract type AbstractElement{dim,M} end

""" Abstract supertype for all elements.

An `AbstractElement` object facilitates the process of evaluating:

    - The internal forces vector and its tangent matrices.
    - The inertial forces vector and its tangent matrices.
    - The mechanical stresses and loads.


**Common methods:**

* [`create_element`](@ref)
* [`coordinates`](@ref)
* [`dofs`](@ref)
* [`local_dofs`](@ref)
* [`index`](@ref)
* [`label`](@ref)
* [`set_label!`](@ref)
* [`nodes`](@ref)
* [`material`](@ref)

* [`internal_forces`](@ref)
* [`inertial_forces`](@ref)

* [`stress`](@ref)
* [`strain`](@ref)

**Common fields:**
* nodes
* material
* label

**Hard contracts:**

* [`create_element`](@ref)  - creates a new element given an empty element defined in the input
* [`local_dofs`](@ref)      - defines the local dofs of the element.

For static cases the following methods are required:

* [`internal_forces`](@ref)     - function that returns the internal force vector.
* [`internal_tangents`](@ref)   - function that returns the internal stiffness matrix.

For dynamic cases the following methods are required:
* [`inertial_forces`](@ref)     - function that returns the inertial force vector.
* [`inertial_tangents`](@ref)   - function that returns the inertial stiffness matrices.

**Soft contracts:**
The following methods can be implemented to provide additional functionality:

* [`stress`](@ref)                - function that returns the element inertial stress.
* [`strain`](@ref)                  - function that returns the element stress.

"""

"Creates a new element given a pre-element (without material assigned) defined in the input"
function create_element(e::AbstractElement, material::AbstractMaterial, args...; kwargs...) end

"Returns element coordinates."
coordinates(e::AbstractElement) = row_vector(coordinates.(nodes(e)))

dofs(e::AbstractElement) = row_vector(dofs.(nodes(e)))
dofs(ve::Vector{<:AbstractElement}) = unique(row_vector(dofs.(ve)))

"Returns local dofs of an element. This dofs are essential for the assemble process."
function local_dofs(e::AbstractElement) end

index(e::AbstractElement) = e.index

"Returns element label."
label(e::AbstractElement) = e.label[]

"Sets element label."
set_label!(e::AbstractElement, label::L) where {L<:Union{String,Symbol}} = e.label[] = Symbol(label)

"Returns element nodes."
nodes(e::AbstractElement) = e.nodes

"Returns element material."
material(e::AbstractElement) = e.material

"Returns the internal force vector of the element."
function internal_forces(e::AbstractElement, args...; kwargs...) end

"Returns the inertial force vector of the element."
function inertial_forces(e::AbstractElement, args...; kwargs...) end

"Returns the element stresses"
function stress(e::AbstractElement, args...; kwargs...) end

"Returns the element strain"
function strain(e::AbstractElement, args...; kwargs...) end


include("./../elements/Truss.jl")

end # module


