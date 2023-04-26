module Utils

using LinearAlgebra: Diagonal

export ScalarWrapper, label, unwrap, eye, row_vector, Point, @debugtime, voigt

#================================#
# Generic functions to overload  #
#================================#
"Return the index of an object "
function index end

"Return the dofs of an object "
function dofs end

"Return the label of an object."
function label end

"Unwraps the object fields."
function _unwrap end

"Scalar mutable struct to avoid making mutable larger structs."
mutable struct ScalarWrapper{T}
    x::T
end

unwrap(s::ScalarWrapper) = s.x
@inline Base.getindex(s::ScalarWrapper) = s.x
@inline Base.setindex!(s::ScalarWrapper, v) = s.x = v

#==================#
# Utils functions  #
#==================#

"Return an eye matrix of size m and type T."
eye(m::Integer, T=Bool) = Diagonal(ones(T, m))

"Transforms a vector of vectors into a 1D row vector."
row_vector(v::Vector{<:AbstractVector{T}}) where {T} = reduce(vcat, v)

"Return the tensor `𝕋` in Voigt notation."
function voigt(𝕋::AbstractMatrix, α::Real=1)
    return [𝕋[1, 1], 𝕋[2, 2], 𝕋[3, 3], α * 𝕋[2, 3], α * 𝕋[1, 3], α * 𝕋[1, 2]]
end

"Execute an expression returning the result and printing the elapsed time inside a `@debug` statement."
macro debugtime(msg, expr)
    quote
        local t0 = time_ns()
        local result = $(esc(expr))
        local t1 = time_ns()
        local elapsed = (t1 - t0) * 1e-9
        local m = $(esc(msg))
        @debug "$m evaluated in $(round(elapsed, digits=8)) sec"
        result
    end
end

end # module
