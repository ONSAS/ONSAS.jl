#################
# Util features #
#################
module Utils

using AutoHashEquals: @auto_hash_equals
using LinearAlgebra: Diagonal

export solve
export dimension, label, set_label!
export ScalarWrapper
export eye, row_vector
export displacements, internal_forces, internal_tangents,
    inertial_forces, inertial_tangents, external_forces, external_tangents

#================================#
# Generic functions to overload  #
#================================#
"Returns the object dimension"
function dimension end

"Returns the object displacements"
function displacements end

"Returns the internal forces of an object"
function internal_forces end

"Returns the inertial forces tangent of an object"
function internal_tangents end

"Returns the inertial forces of an object"
function inertial_forces end

"Returns the inertial forces tangent of an object"
function inertial_tangents end

"Returns the inertial forces of an object"
function external_forces end

"Returns the external forces tangent of an object"
function external_tangents end

"Returns the label of an object."
function label end

"Sets the label of an object."
function set_label! end

"Empty function to solve a problem"
function solve end

#==========#
# Wrappers #
#==========#

"Scalar mutable struct to avoid making mutable larger structs"
mutable struct ScalarWrapper{T}
    x::T
end

@inline Base.getindex(s::ScalarWrapper) = s.x
@inline Base.setindex!(s::ScalarWrapper, v) = s.x = v
Base.copy(s::ScalarWrapper{T}) where {T} = ScalarWrapper{T}(copy(s.x))

#==========================#
# LinearAlgebra functions  #
#==========================#

"Returns an eye matrix of size m and type T."
eye(m::Integer, T=Bool) = Diagonal(ones(T, m))

"Transforms a vector of vectors into a 1D row vector."
row_vector(v::Vector{<:AbstractVector{T}}) where {T} = reduce(vcat, v)

end # module