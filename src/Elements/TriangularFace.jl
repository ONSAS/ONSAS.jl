using ..Elements
using ..Utils
using LinearAlgebra
using Reexport

@reexport import ..Elements: area, create_entity, normal_direction

export TriangularFace, normal_direction

"""
A `TriangularFace` represents an element composed by three `Node`s.
### Fields:
- `nodes`    -- stores triangle nodes.
- `label` -- stores the triangle label.
"""
struct TriangularFace{dim,T<:Real,N<:AbstractNode{dim,T}} <: AbstractFace{dim,T}
    nodes::SVector{3,N}
    label::Label
    function TriangularFace(nodes::SVector{3,N},
                            label::Label=NO_LABEL) where
             {dim,T<:Real,N<:AbstractNode{dim,T}}
        @assert 2 ≤ dim ≤ 3 "TriangularFace is only defined for 2 < dim ≤ 3"
        return new{dim,T,N}(nodes, Symbol(label))
    end
end

"Constructor for a `TriangularFace` element considering the nodes `n₁` `n₂` and `n₃`."
function TriangularFace(n₁::N, n₂::N, n₃::N,
                        label::Label=NO_LABEL) where
         {dim,T<:Real,N<:AbstractNode{dim,T}}
    return TriangularFace(SVector(n₁, n₂, n₃), label)
end

"Constructor for a `TriangularFace` element without nodes and a `label`. This function is used to create meshes via GMSH."
function TriangularFace(label::Label=NO_LABEL)
    return TriangularFace(SVector(Node(0, 0), Node(0, 0), Node(0, 0)), label)
end

"Return the area vector with direction and modulus of a `TriangularFace` element `tf`."
function _area_vec(tf::TriangularFace)
    return cross(coordinates(tf)[2] - coordinates(tf)[1], coordinates(tf)[3] - coordinates(tf)[1]) /
           2
end

"Return the area of a `TriangularFace` element `tf`."
function area(tf::TriangularFace)
    A = norm(_area_vec(tf))
    iszero(A) && error("Area of TriangularFace is zero. Check that nodes are not aligned.")
    return A
end

"Return a `TriangularFace` given an empty `TriangularFace` `tf` and a `Vector` of `Node`s `vn`."
function create_entity(tf::TriangularFace, vn::AbstractVector{<:AbstractNode})
    return TriangularFace(vn[1], vn[2], vn[3], label(tf))
end

"Return the normal direction `n` of a `TriangularFace` element `tf`."
function normal_direction(tf::TriangularFace)
    Atf = _area_vec(tf)
    return Atf / norm(Atf)
end
