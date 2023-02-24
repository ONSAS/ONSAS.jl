using ..Elements: AbstractFace, AbstractNode
using ..Utils: eye
using LinearAlgebra: cross, norm

import ..Elements: area, normal_direction

export TriangularFace, normal_direction

"""
A `TriangularFace` represents an element composed by three `Node`s.
### Fields:
- `nodes`    -- stores triangle nodes.
- `label` -- stores the triangle label.
"""
struct TriangularFace{dim,T<:Real,N<:AbstractNode{dim,T}} <: AbstractFace{dim,T}
    nodes::SVector{3,N}
    label::Symbol
    function TriangularFace(nodes::SVector{3,N}, label=:no_labelled_face) where
    {dim,T<:Real,N<:AbstractNode{dim,T}}
        @assert 2 < dim ≤ 3 "TriangularFace is only defined for 2 < dim ≤ 3"
        new{dim,T,N}(nodes, Symbol(label))
    end
end

"Constructor for a `TriangularFace` element considering the nodes `n₁` `n₂` and `n₃`."
function TriangularFace(n₁::N, n₂::N, n₃::N, label::L=:no_labelled_face) where
{dim,T<:Real,N<:AbstractNode{dim,T},L<:Union{String,Symbol}}
    TriangularFace(SVector(n₁, n₂, n₃), Symbol(label))
end

"Returns the area vector with direction and modulus of a `TriangularFace` element `tf`."
_area_vec(tf::TriangularFace) = cross(coordinates(tf)[2] - coordinates(tf)[1], coordinates(tf)[3] - coordinates(tf)[1]) / 2

"Returns the area of a `TriangularFace` element `tf`."
function area(tf::TriangularFace)
    A = norm(_area_vec(tf))
    iszero(A) && error("Area of TriangularFace is zero. Check that nodes are not aligned.")
    return A
end

"Returns the normal direction `n` of a `TriangularFace` element `tf`."
function normal_direction(tf::TriangularFace)
    Atf = _area_vec(tf)
    return Atf / norm(Atf)
end
