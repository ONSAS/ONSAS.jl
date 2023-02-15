using SparseArrays: sparse, SparseMatrixCSC
using ..Elements: Dof
using ..StructuralModel: AbstractStructure, num_dofs

export Assembler

"Assembler struct to store column indexes, row indexes and values to be inserted in the sparse matrix"
struct Assembler{T}
    I::Vector{Int}
    J::Vector{Int}
    V::Vector{T}
end

function Assembler(N)
    I = Int[]
    J = Int[]
    V = Float64[]
    sizehint!(I, N)
    sizehint!(J, N)
    sizehint!(V, N)

    Assembler(I, J, V)
end

Assembler(s::AbstractStructure) = Assembler(num_dofs(s))

"""Assembles the element matrix `Ke` into `a`."""
function _assemble!(a::Assembler{T}, dofs::AbstractVector{Dof}, Ke::AbstractMatrix{T}) where {T}
    _assemble!(a, index.(dofs), index.(dofs), Ke)
end

""" Assembles the matrix `Ke` into `a` according to the dofs specified by `row_dof_indexes` and `col_dof_indexes`. """
function _assemble!(a::Assembler{T},
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

"Finish the assemble and adds the assembler into the sparse matrix"
function end_assemble!(Kg::AbstractMatrix{T}, a::Assembler{T}) where {T}
    # Fill global K with the assembled values
    @inbounds for i in 1:length(a.I)
        Kg[a.I[i], a.J[i]] += a.V[i]
    end
end


