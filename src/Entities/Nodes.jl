"""
Module defining an interface, `AbstractNode`.

A Node is a point in space with degrees of freedom.
"""
module Nodes

using Reexport
using StaticArrays
using Dictionaries: Dictionary

using ..Utils

@reexport import ..Utils: label, dofs, index

export Dof, Point, AbstractNode, Node, dimension, coordinates, create_node, set_dofs!

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

**Abstract Methods**
* [`coordinates`](@ref)
* [`dimension`](@ref)
* [`dofs`](@ref)

**Abstract fields**
* `x`: coordinates of the node
* `dofs`: mapping from field labels to degrees of freedom

"""
abstract type AbstractNode{dim,T} <: StaticArray{Tuple{dim},T,1} end

"Return node coordinates."
coordinates(n::AbstractNode) = n.x
coordinates(vn::AbstractVector{<:AbstractNode}) = coordinates.(vn)

"Return the node dimension (1D, 2D or 3D)."
dimension(::AbstractNode{dim}) where {dim} = dim

"Return a tuple containing the dimensions of the node."
Base.size(n::AbstractNode) = size(n.x)

"the coordinate at component `i` from the node."
Base.getindex(n::AbstractNode, i::Int) = n.x[i]

"Return node degrees of freedom."
dofs(n::AbstractNode) = n.dofs
dofs(vn::Vector{<:AbstractNode}) = vcat(dofs.(vn)...) # mapreduce(dofs, vcat, vn)
dofs(n::AbstractNode, s::Field) = n.dofs[s]

"Sets a vector of dofs `vd` to the node assigned to the field `s`."
function set_dofs!(n::AbstractNode, s::Field, vd::Vector{Dof})
    unique!(append!(get!(dofs(n), s, Dof[]), vd))
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

"Show the a node."
function Base.show(io::IO, ::MIME"text/plain", n::Node)
    println("• Node at $(n.x) and dofs $(n.dofs)")
end

"1D node constructor with one dim coordinates."
function Node(x₁::T, dofs::Dictionary=Dictionary{Field,Vector{Dof}}()) where {T<:Real}
    Node(Point(x₁), dofs)
end

"2D node constructor with two coordinates."
function Node(x₁::T, x₂::T, dofs::Dictionary=Dictionary{Field,Vector{Dof}}()) where {T<:Real}
    Node(Point(x₁, x₂), dofs)
end

"3D node constructor with coordinates."
function Node(x₁::T, x₂::T, x₃::T,
              dofs::Dictionary=Dictionary{Field,Vector{Dof}}()) where {T<:Real}
    Node(Point(x₁, x₂, x₃), dofs)
end

"Node constructor with a `NTuple`."
function Node(t::NTuple{dim,T},
              dofs::Dictionary=Dictionary{Field,Vector{Dof}}()) where {dim,T<:Real}
    Node(Point(t), dofs)
end

"`Node` constructor with an `AbstractVector` data type."
function Node(v::AbstractVector{T},
              dofs::Dictionary=Dictionary{Field,Vector{Dof}}()) where {T<:Real}
    Node(Point(v...), dofs)
end

"Method to inherit from StaticArrays for a node."
StaticArrays.similar_type(::Type{Node{dim,T}}, ::Type{T}, s::Size{dim}) where {dim,T} = Node{dim,T}
Base.length(::Type{Node{dim,T}}) where {dim,T} = dim

end
