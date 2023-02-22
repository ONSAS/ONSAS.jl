using ..Elements: AbstractFace, AbstractNode
using ..Utils: eye

export TriangularFace

"""
A `TriangularFace` represents an element composed by three `Node`s.
### Fields:
- `n₁`             -- stores first triangle node.
- `n₂`             -- stores second triangle node.
- `n₂`             -- stores third triangle node.
- `label`          -- stores the triangle label.
"""
struct TriangularFace{dim,N<:AbstractNode{dim},T<:Real} <: AbstractFace{dim,T}
    n₁::N
    n₂::N
    n₃::N
    label::Symbol
    function TriangularFace(n₁::N, n₂::N, n₃::N, label=:no_labelled_face) where
    {dim,T<:Real,N<:AbstractNode{dim,T}}
        new{dim,N,T}(n₁, n₂, n₃, Symbol(label))
    end
end

"Returns the nodes for the `TriangularFace` `tf`."
nodes(tf::TriangularFace) = [tf.n₁, tf.n₂, tf.n₃]

