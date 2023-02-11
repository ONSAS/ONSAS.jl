#################
# Util features #
#################
module Utils

using LinearAlgebra: Diagonal

export dimension, label
export ScalarWrapper
export eye, row_vector
export coordinates, nodes, displacements, dofs, index, internal_forces, inertial_forces, external_forces, external_tangents, _unwrap

#================================#
# Generic functions to overload  #
#================================#
"Returns the nodes of an object"
function nodes end

"Returns the coordinates of an object"
function coordinates end

"Returns the object dimension"
function dimension end

"Returns dofs of an object "
function dofs end

"Returns the index of an object "
function index end

"Returns the object displacements"
function displacements end

"Returns the internal forces of an object"
function internal_forces end

"Returns the inertial forces of an object"
function inertial_forces end

"Returns the inertial forces of an object"
function external_forces end

"Returns the label of an object."
function label end

"Unwraps the object fields."
function _unwrap end

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