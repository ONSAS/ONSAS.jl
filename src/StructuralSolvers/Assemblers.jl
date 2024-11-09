module Assemblers

using Dictionaries: Dictionary
using SparseArrays: sparse
using Reexport

using ..Entities
using ..Nodes
using ..Structures

@reexport import ..Entities: elements_cache

export Assembler, assemble!, reset_assembler!, end_assemble, end_assemble!, reset!

"Struct that stores column indexes, row indexes and values for the assemble process."
struct Assembler{T,DT}
    "Column indexes."
    I::Vector{Int}
    "Row indexes."
    J::Vector{Int}
    "Values."
    V::Vector{T}
    "Elements cache."
    cache::DT
end

"Constructor of an `Assembler` with size of the sparse matrix `N` ."
function Assembler(N::Integer, cache=nothing)
    I = Int[]
    J = Int[]
    V = Float64[]
    sizehint!(I, N)
    sizehint!(J, N)
    sizehint!(V, N)

    return Assembler(I, J, V, cache)
end

Assembler(s::AbstractStructure, cache=nothing) = Assembler(num_free_dofs(s), cache)

function elements_cache(a::Assembler, e::AbstractElement)
    a.cache[nameof(typeof(e))]
end

"Assembles the element matrix `Ke` into the `Assembler` struct `a`."
function assemble!(a::Assembler{T}, dofs::AbstractVector{Dof}, Ke::AbstractMatrix{T}) where {T}
    return assemble!(a, dofs, dofs, Ke)
end

"Assembles the matrix `Ke` into `Assembler` `a` according to the dofs specified by `row_dof_indexes` and `col_dof_indexes`."
function assemble!(a::Assembler{T},
                   row_dof_indexes::AbstractVector{Int},
                   col_dof_indexes::AbstractVector{Int},
                   Ke::AbstractMatrix{T}) where {T}
    nrows = length(row_dof_indexes)
    ncols = length(col_dof_indexes)
    @assert size(Ke) == (nrows, ncols) "The size of the element matrix Ke does not match the number of dofs."

    append!(a.V, Ke)
    @inbounds for i in 1:ncols
        append!(a.I, row_dof_indexes)
        for _ in 1:nrows
            push!(a.J, col_dof_indexes[i])
        end
    end
end

"Empties the `Assembler` object `a`. This is useful to reuse the same assembler."
function reset!(a::Assembler{T}) where {T}
    N = length(a.I)
    empty!(a.I)
    empty!(a.J)
    empty!(a.V)
    sizehint!(a.I, N)
    sizehint!(a.J, N)
    return sizehint!(a.V, N)
end

"Return an assembled `AbstractSparseMatrix` from the `Assembler` object `a`."
end_assemble(a::Assembler{T}) where {T} = sparse(a.I, a.J, a.V)

"Inserts an the `Assembler` object `a` into the tangent system matrix `K_sys`."
function end_assemble!(K_sys::AbstractMatrix{T}, a::Assembler{T}) where {T}
    for (index_v, v) in enumerate(a.V)
        K_sys[a.I[index_v], a.J[index_v]] += v
    end
end

end # module
