#################
# Util features #
#################
module Utils

using AutoHashEquals: @auto_hash_equals
using LinearAlgebra: Diagonal

export solve
export dimension, label, elements, nodes, set_label!, index, set_index!
export DofIndex, ElementIndex, NodeIndex, ScalarWrapper
export eye, row_vector

####################################
# Empty functions to be overloaded #
####################################
"Returns the object dimension"
function dimension end

"Returns the degrees of freedom of an object"
function dofs end

"Returns the elements of an object."
function elements end

"Returns the label of an object."
function label end

"Returns the identification number of an object"
function index end

"Returns the nodes of an object."
function nodes end

"Sets the label of an object."
function set_label! end

"Empty function to solve a problem"
function solve end

"Sets an index to an object"
function set_index! end

#############
# Indexes  #
#############

""" Abstract supertype for all indexes.

The following methods are provided by the interface:


**Common methods:**

* [`Base.getindex`](@ref)
* [`Base.setindex!`](@ref)
* [`Base.isequal`](@ref)

"""

abstract type AbstractIndex{I} end

@inline Base.getindex(v::AbstractVector, i::AbstractIndex) = v[i[]]
@inline Base.getindex(i::AbstractIndex) = i.id
@inline Base.setindex!(i::AbstractIndex, id) = i.id = id # callable with i[] = id
@inline Base.isequal(i₁::AbstractIndex, i₂::AbstractIndex) = i₁[] == i₂[]
@inline Base.:(==)(i₁::AbstractIndex, i₂::AbstractIndex) = isequal(i₁, i₂)

"""
`Dof` identification number.
### Fields:
`id` -- integer number. 
"""
@auto_hash_equals mutable struct DofIndex{I<:Integer} <: AbstractIndex{I}
    id::I
end

"""
`Element` identification number.
### Fields:
`id` -- integer number. 
"""
@auto_hash_equals mutable struct ElementIndex{I<:Integer} <: AbstractIndex{I}
    id::I
end

"""
`Node`` identification number.
### Fields:
`id` -- index number. 
"""
@auto_hash_equals mutable struct NodeIndex{I<:Integer} <: AbstractIndex{I}
    id::I
end

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