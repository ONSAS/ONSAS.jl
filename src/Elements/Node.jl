using Dictionaries: Dictionary
using ..Elements: AbstractNode, Dof
export Node
using ..Utils

import StaticArrays

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
        return new{dim,T}(x, dofs)
    end
end

"1D `Node` constructor with coordinates `x₁`."
function Node(x₁::T, dofs::Dictionary=Dictionary{Field,Vector{Dof}}()) where {T<:Real}
    return Node(Point(x₁), dofs)
end

"2D `Node` constructor with coordinates `x₁` and `x₂`."
function Node(x₁::T, x₂::T, dofs::Dictionary=Dictionary{Field,Vector{Dof}}()) where {T<:Real}
    return Node(Point((x₁, x₂)), dofs)
end

"3D `Node` constructor with coordinates `x₁`, `x₂`  and `x₃` and an optional `Dof` dictionary `dofs`.."
function Node(x₁::T, x₂::T, x₃::T,
              dofs::Dictionary=Dictionary{Field,Vector{Dof}}()) where {T<:Real}
    return Node(Point(x₁, x₂, x₃), dofs)
end

"`Node` constructor with a `NTuple` `t` and an optional `Dof` dictionary `dofs`."
function Node(t::NTuple{dim,T},
              dofs::Dictionary=Dictionary{Field,Vector{Dof}}()) where {dim,T<:Real}
    return Node(Point(t), dofs)
end

"`Node` constructor with an `AbstractVector` `v`."
function Node(v::AbstractVector{T},
              dofs::Dictionary=Dictionary{Field,Vector{Dof}}()) where {T<:Real}
    return Node(Point(v...), dofs)
end

StaticArrays.similar_type(::Type{Node{dim,T}}, ::Type{T}, s::Size{dim}) where {dim,T} = Node{dim,T}

Base.length(::Type{Node{dim,T}}) where {dim,T} = dim
