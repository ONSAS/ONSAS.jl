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

####################################
# Empty functions to be overloaded #
####################################
"Returns the object dimension"
function dimension end

"Returns the label of an object."
function label end

"Sets the label of an object."
function set_label! end

"Empty function to solve a problem"
function solve end

#############
# Wrappers  #
#############

"Scalar mutable struct to avoid making mutable larger structs"
mutable struct ScalarWrapper{T}
    x::T
end

@inline Base.getindex(s::ScalarWrapper) = s.x
@inline Base.setindex!(s::ScalarWrapper, v) = s.x = v
Base.copy(s::ScalarWrapper{T}) where {T} = ScalarWrapper{T}(copy(s.x))

###################################
# Useful LinearAlgebra functions  #
###################################

"Returns an eye matrix of size m and type T."
eye(m::Integer, T=Bool) = Diagonal(ones(T, m))

"Transforms a vector of vectors into a 1D row vector."
row_vector(v::Vector{<:AbstractVector{T}}) where {T} = reduce(vcat, v)

end # module