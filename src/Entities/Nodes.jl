"""
Module defining an interface, `AbstractNode`.

A Node is a point in space with degrees of freedom.
"""
module Nodes

using Reexport
using StaticArrays
using Dictionaries: Dictionary

using ..Utils

@reexport import ..Utils: label, apply!, dofs, index

export Dof, Point, AbstractNode, Node, dimension, coordinates, create_node

"""
Scalar degree of freedom of the structure.
"""
const Dof = Int

"the dof index of the `Dof` `d` "
index(d::Dof) = d

"""
Construct a point given its coordinates.

It is a type alias for a statically sized array of dimension `dim` and element type `T`.
Use either as a vararg function, `Point(1, 2, 3)`, or by splatting `Point(x...)` if `x` is an (abstract) vector.
"""
const Point{dim,T} = SVector{dim,T} where {dim,T<:Real}

"""
An `AbstractNode` object is a point in space with degrees of freedom.

### Methods

The following methods are provided by the interface:

* [`coordinates`](@ref)
* [`dimension`](@ref)
* [`dofs`](@ref)
"""
abstract type AbstractNode{dim,T} <: StaticArray{Tuple{dim},T,1} end

"the `AbstractNode` `n` coordinates."
coordinates(n::AbstractNode) = n.x

"each `AbstractNode` coordinates in a `Vector` of `Node`s vn."
coordinates(vn::AbstractVector{<:AbstractNode}) = coordinates.(vn)

"the `AbstractNode` `n` dimension (1D, 2D or 3D)."
dimension(::AbstractNode{dim}) where {dim} = dim

"a tuple containing the dimensions of the `AbstractNode` `n`."
Base.size(n::AbstractNode) = size(n.x)

"the coordinate at component `i` from the `AbstractNode` `n`."
Base.getindex(n::AbstractNode, i::Int) = n.x[i]

"`AbstractNode` `n` degrees of freedom."
dofs(n::AbstractNode) = n.dofs

"degrees of freedom for each `AbstractNode` in a vector of nodes `vn`."
dofs(vn::Vector{<:AbstractNode}) = vcat(dofs.(vn)...)

"`AbstractNode` `n` degrees of freedom with symbol `s`."
dofs(n::AbstractNode, s::Field) = n.dofs[s]

"Sets a `Vector`s of dofs `vd` to the `AbstractNode` `n` assigned to the field `s`."
function apply!(n::AbstractNode, s::Field, vd::Vector{Dof})
    if s ∉ keys(dofs(n))
        insert!(dofs(n), s, vd)
    else
        [push!(dofs(n)[s], d) for d in vd if d ∉ dofs(n)[s]]
    end
    dofs(n)
end

"""
A `Node` is a point in space.
The coordinates of the node are stored using a static array.
"""
struct Node{dim,T} <: AbstractNode{dim,T}
    "Coordinates of the node."
    x::Point{dim,T}
    "Mapping from field labels to degrees of freedom."
    dofs::Dictionary{Field,Vector{Dof}}
    function Node(x::Point{dim,T}, dofs::Dictionary{Field,Vector{Dof}}) where {dim,T<:Real}
        @assert dim ≤ 3 "Only 1D, 2D or 3D nodes are supported."
        new{dim,T}(x, dofs)
    end
end

"Show the `Node` `n`."# add mimeTexPlain
function Base.show(io::IO, ::MIME"text/plain", n::Node)
    println("• Node at $(n.x) and dofs $(n.dofs)")
end

"1D `Node` constructor with coordinates `x₁`."
function Node(x₁::T, dofs::Dictionary=Dictionary{Field,Vector{Dof}}()) where {T<:Real}
    Node(Point(x₁), dofs)
end

"2D `Node` constructor with coordinates `x₁` and `x₂`."
function Node(x₁::T, x₂::T, dofs::Dictionary=Dictionary{Field,Vector{Dof}}()) where {T<:Real}
    Node(Point((x₁, x₂)), dofs)
end

"3D `Node` constructor with coordinates `x₁`, `x₂`  and `x₃` and an optional `Dof` dictionary `dofs`.."
function Node(x₁::T, x₂::T, x₃::T,
              dofs::Dictionary=Dictionary{Field,Vector{Dof}}()) where {T<:Real}
    Node(Point(x₁, x₂, x₃), dofs)
end

"`Node` constructor with a `NTuple` `t` and an optional `Dof` dictionary `dofs`."
function Node(t::NTuple{dim,T},
              dofs::Dictionary=Dictionary{Field,Vector{Dof}}()) where {dim,T<:Real}
    Node(Point(t), dofs)
end

"`Node` constructor with an `AbstractVector` `v`."
function Node(v::AbstractVector{T},
              dofs::Dictionary=Dictionary{Field,Vector{Dof}}()) where {T<:Real}
    Node(Point(v...), dofs)
end

"Method to inherit from StaticArrays for a `Node`."
StaticArrays.similar_type(::Type{Node{dim,T}}, ::Type{T}, s::Size{dim}) where {dim,T} = Node{dim,T}
Base.length(::Type{Node{dim,T}}) where {dim,T} = dim

end
