module Utils

using LinearAlgebra

export label, eye, row_vector, @debugtime, voigt, Label, NO_LABEL, Density, Field, index,
       fill_symmetric_matrix!, INDEXES_TO_VOIGT

#================================#
# Generic functions to overload  #
#================================#
"Return the index of an object "
function index end

"Return the dofs of an object "
function dofs end

"Return the label of an object."
function label end

"Apply one object to the other."
function apply! end

#==================#
# Utils functions  #
#==================#

"Return an eye matrix of size m and type T."
eye(m::Integer, T=Bool) = Diagonal(ones(T, m))

"Transforms a vector of vectors into a 1D row vector."
row_vector(v::Vector{<:AbstractVector{T}}) where {T} = reduce(vcat, v)

function fill_symmetric_matrix!(S::Symmetric{T},
                                S₁₁::T, S₂₂::T, S₃₃::T, S₂₃::T, S₁₃::T, S₁₂::T) where {T}
    A = parent(S)
    S[1, 1] = S₁₁
    S[2, 2] = S₂₂
    S[3, 3] = S₃₃
    A[2, 3] = S₂₃
    A[1, 3] = S₁₃
    A[1, 2] = S₁₂
end

"Indexes to transform form Tensors.jl to Voigt nomenclature."
const CONSTANT_INDICES_TO_VOIGT = [[1, 1], [2, 2], [3, 3]]
const SCALED_INDICES_TO_VOIGT = [[2, 3], [1, 3], [1, 2]]

#TODO: Replace indexes
"Return the tensor `𝕋` in Voigt notation."
function voigt(𝕋::AbstractMatrix, α::Real=1)
    @assert size(𝕋) == (3, 3) "Unexpected tensor dimensions"
    v = [𝕋[idx...] for idx in CONSTANT_INDICES_TO_VOIGT]
    for idx in SCALED_INDICES_TO_VOIGT
        push!(v, α * 𝕋[idx]...)
    end
    v
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

#==================#
# Type aliases     #
#==================#

"Used to assign labels to geometric or physical entities."
const Label = Union{String,Symbol}

"Label to design an entity without assigned label."
const NO_LABEL = :no_label

"Physical paramater defining density of a material (`nothing` reserved for static cases)."
const Density = Union{Float64,Nothing}

"Type alias used for field labels as degree-of-freedom keys such as `:T` or `:θ`."
const Field = Symbol

end # module
