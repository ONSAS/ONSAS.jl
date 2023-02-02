######################
# Elements interface #
######################

"""
Module defining the elements implemented.
"""
module Elements

using ..Materials: AbstractMaterial
using ..CrossSections: AbstractCrossSection, area
using ..BoundaryConditions: AbstractBoundaryCondition
using ..Utils: ScalarWrapper
using Reexport: @reexport
using StaticArrays: SVector, SMatrix

@reexport import ..Utils: Index, index, dofs, dimension, nodes, label, set_label!, set_index!

export Dof, symbol, is_fixed, fix!
export AbstractNode, Node, boundary_conditions, coordinates, coordinates_eltype, set_dof_index!
export AbstractElement, num_nodes, dofs_per_node, geometry, material, material_model
export Truss, internal_force, stiffness_matrix


const _DEFAULT_LABEL = :no_labelled_element
const _DEFAULT_INDEX_INT = 0
const _DEFAULT_INDEX = Index(_DEFAULT_INDEX_INT)

""" Degree of freedom struct.
This is a scalar degree of freedom of the structure.
### Fields:
- `symbol`  -- degree of freedom symbol.
- `index`   -- degree of freedom identification number. 
- `is_free` -- boolean indicating if the dof is free(`false`) or fixed(`true`).
"""
struct Dof
    symbol::Symbol
    index::Index
    is_fixed::ScalarWrapper{Bool}
end

Dof(symbol::Symbol, index::Integer=_DEFAULT_INDEX_INT, is_fixed::Bool=false) = Dof(symbol, Index(index), ScalarWrapper(is_fixed))
Base.:(==)(d1::Dof, d2::Dof) = d1.symbol == d2.symbol && d1.index == d2.index


"Returns the degree of freedom identification number."
index(d::Dof) = d.index

"Sets a new index to the degree of freedom."
set_index!(d::Dof, i::Int) = d.index[] = i
set_index!(d::Dof, idx::Index) = d.index = idx

"Returns the degree of freedom symbol."
symbol(d::Dof) = d.symbol

"Returns `true` if the degree of freedom is free"
is_fixed(d::Dof) = d.is_fixed[]

"Sets the degree of freedom as fixed"
fix!(d::Dof) = d.is_fixed[] = true

# =================
# Abstract Node
# =================

abstract type AbstractNode{dim,T} end

""" Abstract supertype for all nodes.

An `AbstractNode` object is a point in space.

**Common methods:**

* [`boundary_conditions`](@ref)
* [`coordinates`](@ref)
* [`coordinates_eltype`](@ref)
* [`index`](@ref)
* [`dimension`](@ref)
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

"Sets node's index."
set_index!(n::AbstractNode, id::Int) = n.index[] = id

"Sets node dofs indexes."
function set_dof_index!(n::AbstractNode, indexes::I) where {I<:Union{Vector{<:Integer},AbstractRange{<:Integer}}}
    ndofs = dofs(n)
    length(ndofs) != length(indexes) && throw(ArgumentError("Indexes and dofs length mismatch"))
    [set_index!(dof, indexes[i]) for (i, dof) in enumerate(ndofs)]
    return ndofs
end

"Fixes a node dofs"
fix!(n::AbstractNode) = [fix!(dof) for dof in dofs(n)]


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
"""
mutable struct Node{dim,T} <: AbstractNode{dim,T}
    x::AbstractArray{T}
    dofs::Vector{<:Dof}
    index::Index
    bc::Vector{<:AbstractBoundaryCondition}
    function Node(x::AbstractArray{T}, dofs::Vector{<:Dof}, index::Index, bc::Vector{<:AbstractBoundaryCondition}) where {T}
        new{length(x),T}(x, dofs, index, bc)
    end
end

Node(x::NTuple{dim,T}, index::Integer, bc::Vector{<:AbstractBoundaryCondition}) where {dim,T} = Node(SVector(x), _dim_to_nodal_dofs(length(x)), Index(index), bc)
# Empty bc constructors
empty_bc = Vector{AbstractBoundaryCondition}()
Node(x::NTuple{dim,T}, index::Index) where {dim,T} = Node(SVector(x), _dim_to_nodal_dofs(length(x)), index, empty_bc)
Node(x::NTuple{dim,T}, index::Integer=_DEFAULT_INDEX_INT) where {dim,T} = Node(SVector(x), _dim_to_nodal_dofs(length(x)), Index(index), empty_bc)
Node(x::AbstractArray{T}, index::Index) where {T} = Node(x, _dim_to_nodal_dofs(length(x)), index, empty_bc)
Node(x::AbstractArray{T}, index::Integer=_DEFAULT_INDEX_INT) where {T} = Node(x, _dim_to_nodal_dofs(length(x)), Index(index), empty_bc)


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
* [`label`](@ref)
* [`nodes`](@ref)
* [`set_nodes!`](@ref)
* [`set_label!`](@ref)

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

coordinates(e::AbstractElement) = row_vector(coordinates.(nodes(e)))

dofs(e::AbstractElement) = row_vector(dofs.(nodes(e)))

"Returns element number of dofs per element"
dofs_per_node(e::AbstractElement) = length(dofs(first(nodes(e))))

dimension(::AbstractElement{dim}) where {dim} = dim

"Returns the geometrical properties of the element"
geometry(e::AbstractElement) = e.geometry

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
function internal_force(e::AbstractElement, args...; kwargs...) end

"Returns the internal stiffness matrix of the element."
function stiffness_matrix(e::AbstractElement, args...; kwargs...) end

"Returns the inertial force vector of the element."
function inertial_force(e::AbstractElement, args...; kwargs...) end

"Returns the inertial tangent matrices of the element."
function mass_matrices(e::AbstractElement, args...; kwargs...) end

include("./../elements/Truss.jl")

end # module


