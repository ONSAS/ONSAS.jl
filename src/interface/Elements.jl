"""
Module defining the elements implemented.
"""
module Elements

using AutoHashEquals: @auto_hash_equals
using Reexport: @reexport
@reexport using Dictionaries
using StaticArrays: SVector
using ..Utils: row_vector

@reexport using ..CrossSections
@reexport import ..Utils: dimension, label
@reexport import ..Utils: internal_forces, inertial_forces
@reexport import ..Utils: coordinates, dofs, nodes
import Dictionaries: index

export Dof, add_dofs!
export AbstractNode, Node, coordinates
export AbstractElement, cross_section, coordinates, local_dof_symbol, local_dofs

# ========================
# Degree of freedom (Dof)
# ========================

""" Degree of freedom struct.
This is a scalar degree of freedom of the structure.
### Fields:
- `index`   -- degree of freedom identification number. 
"""
@auto_hash_equals struct Dof
    index::Int
end

index(d::Dof) = d.index
index(vd::Vector{Dof}) = index.(vd)
Base.maximum(vd::Vector{Dof}) = maximum(index.(vd))

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
dofs(n::AbstractNode, s::Symbol) = n.dofs[s]
dofs(n::Vector{<:AbstractNode}) = vcat(dofs.(n)...)

"Adds a vectors of `Dof` with a symbol"
function add_dofs!(n::AbstractNode, s::Symbol, dofs_to_add::Vector{Dof})
    if s ∉ keys(dofs(n))
        insert!(dofs(n), s, dofs_to_add)
    else
        [push!(dofs(n)[s], d) for d in dofs_to_add if d ∉ dofs(n)[s]]
    end
    return dofs(n)
end


"""
A `Node` is a point in space.
### Fields:
- `x`     -- stores the coordinates.
- `dofs`  -- stores the node degrees of freedom, maps symbol to dofs.
"""
struct Node{dim,T} <: AbstractNode{dim,T}
    x::AbstractArray{T}
    dofs::Dictionary{Symbol,Vector{Dof}}
    function Node(
        x::AbstractArray{T}, dofs::Dictionary{Symbol,Vector{Dof}}=Dictionary{Symbol,Vector{Dof}}()) where {T<:Real}
        dim = length(x)
        @assert dim ≤ 3 "Only 1D,2D or 3D nodes are supported"
        new{dim,T}(x, dofs)
    end
end

Node(x₁::T) where {T<:Real} = Node(SVector(x₁))
Node(x₁::T, x₂::T) where {T<:Real} = Node(SVector((x₁, x₂)))
Node(x₁::T, x₂::T, x₃::T) where {T<:Real} = Node(SVector((x₁, x₂, x₃)))
Node(x::NTuple{dim,T}) where {dim,T<:Real} = Node(SVector(x))
Node(x::Vector{T}) where {T<:Real} = Node(SVector(x...))

# =================
# Abstract Element
# =================

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

* [`local_dof_symbol`](@ref)      - defines the local dofs of the element.

For static cases the following methods are required:

* [`inertial_forces`](@ref) - function that returns the internal forces vector and its respective tangent matrices.

For dynamic cases the following methods are required:
* [`inertial_forces`](@ref) - function that returns the inertial forces vector and its respective tangent matrices.
"""

"Returns element coordinates."
coordinates(e::AbstractElement) = coordinates.(nodes(e))

"Returns the geometrical properties of the element"
cross_section(e::AbstractElement) = e.cross_section

function dofs(e::AbstractElement)

    element_dofs = Dictionary{Symbol,Vector{Dof}}()
    nodes_dofs_vec = dofs.(nodes(e))

    for node_dofs in nodes_dofs_vec
        for dof_symbol in keys(node_dofs)
            dofs_to_add = node_dofs[dof_symbol]
            if dof_symbol ∉ keys(element_dofs)
                insert!(element_dofs, dof_symbol, dofs_to_add)
            else
                [push!(element_dofs[dof_symbol], d) for d in dofs_to_add if d ∉ element_dofs[dof_symbol]]
            end
        end
    end
    return element_dofs
end

dofs(ve::Vector{<:AbstractElement}) = unique(row_vector(dofs.(ve)))

"Returns local dofs symbols of the element (for only displacements `:u` is used) in a vector. 
This dofs are essential for the assemble process."
function local_dof_symbol(e::AbstractElement) end

"Returns local dofs given a vector of local dof symobls. This method extracts all node dofs with the same symbol
as local_dof_symbol"
function local_dofs(e::AbstractElement)
    local_dof_symbols = local_dof_symbol(e)
    local_dofs = Vector{Dof}()
    element_dofs = dofs(e)
    for dof_symbol in local_dof_symbols
        if dof_symbol ∉ keys(element_dofs)
            error("Element $(e.label) does not have dofs with symbol $(dof_symbol)")
        else
            push!(local_dofs, element_dofs[dof_symbol]...)
        end
    end
    return local_dofs
end


"Returns element label."
label(e::AbstractElement) = e.label

"Returns element nodes."
nodes(e::AbstractElement) = e.nodes

"Returns the internal force vector of the element."
function internal_forces(e::AbstractElement, args...; kwargs...) end

"Returns the inertial force vector of the element."
function inertial_forces(e::AbstractElement, args...; kwargs...) end

include("../elements/Truss.jl")

end # module


