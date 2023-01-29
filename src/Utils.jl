#################
# Util features #
#################

module Utils

using AutoHashEquals: @auto_hash_equals
export Index, ScalarWrapper, index, set_index!, label, solve

"Empty function to extract the label of an object."
function label end

"Empty function to solve a problem"
function solve end

"Returns the identification number of an object"
function index end

"Sets an index to an object"
function set_index! end

"Returns the degrees of freedom of an object"
function dofs end

"Scalar mutable struct to avoid using larger mutable structs"
@auto_hash_equals mutable struct Index
    id::Integer
end

@inline Base.getindex(i::Index) = i.id
@inline Base.setindex!(i::Index, id) = i.id = id # callable with i[] = id


"Scalar mutable struct to avoid making mutable larger structs"
mutable struct ScalarWrapper{T}
    x::T
end

@inline Base.getindex(s::ScalarWrapper) = s.x
@inline Base.setindex!(s::ScalarWrapper, v) = s.x = v
Base.copy(s::ScalarWrapper{T}) where {T} = ScalarWrapper{T}(copy(s.x))

end # module