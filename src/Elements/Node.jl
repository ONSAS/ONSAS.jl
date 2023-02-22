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
struct Node{dim,T,X} <: AbstractNode{dim,T}
    x::X
    dofs::Dictionary{Symbol,Vector{Dof}}
    function Node(
        x::X, dofs::Dictionary{Symbol,Vector{Dof}}=Dictionary{Symbol,Vector{Dof}}()) where {T<:Real,X<:AbstractVector{T}}
        dim = length(x)
        @assert dim ≤ 3 "Only 1D,2D or 3D nodes are supported"
        new{dim,T,X}(x, dofs)
    end
end

"1D `Node` constructor with coordinates `x₁`."
Node(x₁::T) where {T<:Real} = Node(SVector(x₁))

"2D `Node` constructor with coordinates `x₁` and `x₂`."
Node(x₁::T, x₂::T) where {T<:Real} = Node(SVector((x₁, x₂)))

"3D `Node` constructor with coordinates `x₁`, `x₂`  and `x₃`."
Node(x₁::T, x₂::T, x₃::T) where {T<:Real} = Node(SVector((x₁, x₂, x₃)))

"`Node` constructor with a `NTuple` `t`."
Node(t::NTuple{dim,T}) where {dim,T<:Real} = Node(SVector(t))

"`Node` constructor with a `Vector` `v`."
Node(v::Vector{T}) where {T<:Real} = Node(SVector(v...))