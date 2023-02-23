using ..Elements: AbstractFace, AbstractNode
using ..Utils: eye
using LinearAlgebra: cross, norm

import ..Elements: area, normal_direction

export TriangularFace, normal_direction

"""
A `TriangularFace` represents an element composed by three `Node`s.
### Fields:
- `n₁`    -- stores first triangle node.
- `n₂`    -- stores second triangle node.
- `n₂`    -- stores third triangle node.
- `label` -- stores the triangle label.
"""
struct TriangularFace{dim,N<:AbstractNode{dim},T<:Real} <: AbstractFace{dim,T}
    nodes::SVector{3,N}
    label::Symbol
    function TriangularFace(n₁::N, n₂::N, n₃::N, label=:no_labelled_face) where
    {dim,T<:Real,N<:AbstractNode{dim,T}}
        @assert dim > 1 "TriangularFace is only defined for dim > 1"
        new{dim,N,T}(SVector(n₁, n₂, n₃), Symbol(label))
    end
end

"Returns the area vector with direction and modulus of a `TriangularFace` element `tf`."
_area_vec(tf::TriangularFace) = cross(coordinates(tf)[2] - coordinates(tf)[1], coordinates(tf)[3] - coordinates(tf)[1]) / 2

"Returns the area of a `TriangularFace` element `tf`."
area(tf::TriangularFace) = norm(_area_vec(tf))

"Returns the normal direction `n` of a `TriangularFace` element `tf`."
function normal_direction(tf::TriangularFace)
    Atf = _area_vec(tf)
    return Atf / norm(Atf)
end

