#################
# Util features #
#################

module Utils

export Index, index, label, solve

"Empty function to extract the label of an object."
function label end

"Empty function to solve a problem"
function solve end

"Returns the identification number of an object"
function index end


"Scalar mutable struct to avoid using larger mutable structs"
mutable struct Index
    id::Integer
end

@inline Base.getindex(i::Index) = i.id
@inline Base.setindex!(i::Index, id) = i.id = id # callable with i[] = id



end # module