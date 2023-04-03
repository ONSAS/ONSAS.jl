#################
# Util features #
#################
module Utils

using LinearAlgebra: Diagonal

export label, _unwrap

export eye, row_vector
#================================#
# Generic functions to overload  #
#================================#
"Returns the index of an object "
function index end

"Returns the dofs of an object "
function dofs end

"Returns the label of an object."
function label end

"Unwraps the object fields."
function _unwrap end

#==================#
# Utils functions  #
#==================#

"Returns an eye matrix of size m and type T."
eye(m::Integer, T=Bool) = Diagonal(ones(T, m))

"Transforms a vector of vectors into a 1D row vector."
row_vector(v::Vector{<:AbstractVector{T}}) where {T} = reduce(vcat, v)

"Returns the Voigt notation of tensor `𝕋`."
_voigt(𝕋::AbstractMatrix, α::Real=1) = [𝕋[1, 1], 𝕋[2, 2], 𝕋[3, 3], α * 𝕋[2, 3], α * 𝕋[1, 3], α * 𝕋[1, 2]]


end # module