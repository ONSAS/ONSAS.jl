######################
# Elements interface #
######################

"""
Module defining the elements implemented.
"""
module Elements

using AutoHashEquals: @auto_hash_equals
using Reexport: @reexport
using StaticArrays: SVector, SMatrix

using ..Materials: AbstractMaterial
using ..CrossSections: AbstractCrossSection, area
using ..BoundaryConditions: AbstractBoundaryCondition
@reexport import ..BoundaryConditions: dofs
using ..InitialConditions: AbstractInitialCondition
using ..Utils: ScalarWrapper
@reexport import ..Utils: dimension, label, set_label!
@reexport import ..Utils: internal_forces, internal_tangents, displacements

export AbstractIndex, ElementIndex, NodeIndex, DofIndex, index, set_index!
export Dof, symbol, is_fixed, fix!
export AbstractNode, Node, boundary_conditions, coordinates, coordinates_eltype, element_type, dofs
export AbstractElement, num_nodes, dofs_per_node, geometry, material, material_model

const _DEFAULT_LABEL = :no_labelled_element

#############
# Indexes  #
#############

""" Abstract supertype for all indexes.

The following methods are provided by the interface:


**Common methods:**

* [`Base.getindex`](@ref)
* [`Base.setindex!`](@ref)
* [`Base.isequal`](@ref)

"""

abstract type AbstractIndex{I} end

@inline Base.getindex(v::AbstractVector, i::AbstractIndex) = v[i[]]
@inline Base.getindex(i::AbstractIndex) = i.id
@inline Base.setindex!(i::AbstractIndex, id) = i.id = id # callable with i[] = id
@inline Base.isequal(i₁::AbstractIndex, i₂::AbstractIndex) = i₁[] == i₂[]
@inline Base.:(==)(i₁::AbstractIndex, i₂::AbstractIndex) = isequal(i₁, i₂)

"""
`Dof` identification number.
### Fields:
`id` -- integer number. 
"""
@auto_hash_equals mutable struct DofIndex{I<:Integer} <: AbstractIndex{I}
    id::I
end

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
const _DEFAULT_DOF_INDEX = DofIndex(_DEFAULT_INDEX_INT)

""" Degree of freedom struct.
This is a scalar degree of freedom of the structure.
### Fields:
- `symbol`  -- degree of freedom symbol.
- `index`   -- degree of freedom identification number. 
- `is_free` -- boolean indicating if the dof is free(`false`) or fixed(`true`).
"""
struct Dof
    symbol::Symbol
    index::DofIndex
    is_fixed::ScalarWrapper{Bool}
end

Dof(symbol::Symbol, index::Integer=_DEFAULT_INDEX_INT, is_fixed::Bool=false) = Dof(symbol, DofIndex(index), ScalarWrapper(is_fixed))
Base.:(==)(d1::Dof, d2::Dof) = d1.symbol == d2.symbol && d1.index == d2.index


"Returns the degree of freedom identification number."
index(d::Dof) = d.index

"Sets a new index to the degree of freedom."
set_index!(d::Dof, i::Int) = d.index[] = i
set_index!(d::Dof, idx::DofIndex) = d.index = idx

"Returns the degree of freedom symbol."
symbol(d::Dof) = d.symbol

"Returns `true` if the degree of freedom is free"
is_fixed(d::Dof) = d.is_fixed[]

"Sets the degree of freedom as fixed"
fix!(d::Dof) = d.is_fixed[] = true

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

function set_index!(n::AbstractNode{dim}, idx::NodeIndex) where {dim}
    ndofs = dofs(n)
    node_dof_indexes = _nodes2dofs(idx[], dim)
    [set_index!(dof, node_dof_indexes[i]) for (i, dof) in enumerate(ndofs)]
    n.index[] = idx[]
    return n
end

set_index!(n::AbstractNode, i::Integer) = set_index!(n, NodeIndex(i))

"Fixes a node dofs"
function fix!(n::AbstractNode)
    [fix!(dof) for dof in dofs(n)]
    return n
end
#TODO: generalize to any field type and dimension (:u, dim)
"Maps dimension of the node to local degrees of freedom."
function _dim_to_nodal_dofs(dim::Int)
    dofs = if dim == 1
        [Dof(:uₓ, 1)]
    elseif dim == 2
        [Dof(:uₓ, 1), Dof(:uⱼ, 2), Dof(:θₖ, 3)]
    elseif dim == 3
        [
            Dof(:uₓ, 1), Dof(:θₓ, 2),
            Dof(:uⱼ, 3), Dof(:θⱼ, 4),
            Dof(:uₖ, 5), Dof(:θₖ, 6)
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

* [`coordinates`](@ref)
* [`dofs`](@ref)
* [`dofs_per_node`](@ref)
* [`dimension`](@ref)
* [`geometry`](@ref)
* [`index`](@ref)
* [`label`](@ref)
* [`nodes`](@ref)
* [`set_nodes!`](@ref)
* [`set_label!`](@ref)
* [`set_index!`](@ref)
* [`element_type`](@ref)

* [`material`](@ref)
* [`material_model`](@ref)
* [`set_material!`](@ref)

* [`internal_forces`](@ref)
* [`stiffness_matrix`](@ref)
* [`inertial_forces`](@ref)
* [`inertial_tangents`](@ref)

* [`strain`](@ref)
* [`stresses`](@ref)


**Common fields:**
* index
* nodes
* material
* geometry

**Hard contracts:**

* [`num_nodes`](@ref)  - defines the number of nodes per element.
* [`local_dofs`](@ref)  - defines the local dofs of the element.

For static cases the following methods are required:

* [`internal_forces`](@ref)  - function that returns the internal force vector.
* [`stiffness_matrix`](@ref) - function that returns the internal stiffness matrix.

For dynamic cases the following methods are required:
* [`inertial_forces`](@ref) - function that returns the inertial force vector.
* [`inertial_tangents`](@ref)  - function that returns the inertial stiffness matrices.

**Soft contracts:**
The following methods can be implemented to provide additional functionality:

* `stresses`                - function that returns the element inertial stress.
* `strain`                  - function that returns the element stress.

"""
Base.push!(e::AbstractElement, bc::AbstractBoundaryCondition) = push!(e.bc, bc)

coordinates(e::AbstractElement) = row_vector(coordinates.(nodes(e)))

dofs(e::AbstractElement) = row_vector(dofs.(nodes(e)))

"Returns element number of dofs per element"
dofs_per_node(e::AbstractElement) = length(dofs(first(nodes(e))))

dimension(::AbstractElement{dim}) where {dim} = dim

"Returns the geometrical properties of the element"
geometry(e::AbstractElement) = e.geometry

index(e::AbstractElement) = e.index

set_index!(e::AbstractElement, i::Int) = set_index(e, ElementIndex(i))

set_index!(e::AbstractElement, idx::ElementIndex) = e.index = idx

nodes(e::AbstractElement) = e.nodes

"Sets nodes to the element."
function set_nodes!(e::AbstractElement, nodes::Vector{<:AbstractNode})
    if length(nodes) == num_nodes(e)
        e.nodes = nodes
    else
        return ArgumentError("The number of nodes must be $(num_nodes(e)).")
    end
end

label(e::AbstractElement) = e.label[]

set_label!(e::AbstractElement, label::String) = set_label!(e, Symbol(label))

set_label!(e::AbstractElement, label::Symbol) = e.label[] = label

"Returns the element material."
material(e::AbstractElement) = e.material

"Returns the material model of the element."
material_model(::AbstractElement{dim,M}) where {dim,M} = M

"Sets a material to the element."
set_material!(e::AbstractElement, m::AbstractMaterial) = e.material = m

"Returns the internal force vector of the element."
function internal_forces(e::AbstractElement, args...; kwargs...) end

"Returns the internal stiffness matrix of the element."
function stiffness_matrix(e::AbstractElement, args...; kwargs...) end

"Returns the inertial force vector of the element."
function inertial_forces(e::AbstractElement, args...; kwargs...) end

"Returns the inertial tangent matrices of the element."
function mass_matrices(e::AbstractElement, args...; kwargs...) end

include("./../elements/Truss.jl")

end # module


