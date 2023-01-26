#################
# Util features #
#################

module Utils

export ScalarWrapper
export label, solve

"Scalar mutable struct to avoid making mutable larger structs"
mutable struct ScalarWrapper{T}
    x::T
end

@inline Base.getindex(s::ScalarWrapper) = s.x
@inline Base.setindex!(s::ScalarWrapper, v) = s.x = v
Base.copy(s::ScalarWrapper{T}) where {T} = ScalarWrapper{T}(copy(s.x))

"Empty function to extract the object label"
function label end

"Empty function to solve a problem"
function solve end


end # module