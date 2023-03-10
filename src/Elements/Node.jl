using Dictionaries: Dictionary
using ..Elements: AbstractNode, Dof
using StaticArrays: SVector
export Node

"""
A `Node` is a point in space. The coordinates of the node are stored using `SVector` data type from static arrays. 
### Fields:
- `x`     -- stores the coordinates.
- `dofs`  -- stores the node degrees of freedom, maps symbol to dofs.
"""
struct Node{dim,T} <: AbstractNode{dim,T}
    x::SVector{dim,T}
    dofs::Dictionary{Symbol,Vector{Dof}}
    function Node(
        x::SVector{dim,T}, dofs::Dictionary{Symbol,Vector{Dof}}) where {dim,T<:Real}
        @assert dim ≤ 3 "Only 1D,2D or 3D nodes are supported"
        new{dim,T}(x, dofs)
    end
end

"1D `Node` constructor with coordinates `x₁`."
Node(x₁::T, dofs::Dictionary=Dictionary{Symbol,Vector{Dof}}()) where {T<:Real} =
    Node(SVector(x₁), dofs)

"2D `Node` constructor with coordinates `x₁` and `x₂`."
Node(x₁::T, x₂::T, dofs::Dictionary=Dictionary{Symbol,Vector{Dof}}()) where {T<:Real} =
    Node(SVector((x₁, x₂)), dofs)


"3D `Node` constructor with coordinates `x₁`, `x₂`  and `x₃` and an optional `Dof` dictionary `dofs`.."
Node(x₁::T, x₂::T, x₃::T, dofs::Dictionary=Dictionary{Symbol,Vector{Dof}}()) where {T<:Real} =
    Node(SVector(x₁, x₂, x₃), dofs)

"`Node` constructor with a `NTuple` `t` and an optional `Dof` dictionary `dofs`."
Node(t::NTuple{dim,T}, dofs::Dictionary=Dictionary{Symbol,Vector{Dof}}()) where {dim,T<:Real} =
    Node(SVector(t), dofs)

"`Node` constructor with an `AbstractVector` `v`."
Node(v::AbstractVector{T}, dofs::Dictionary=Dictionary{Symbol,Vector{Dof}}()) where {T<:Real} =
    Node(SVector(v...), dofs)