using ..StructuralModel: AbstractStructure, num_dofs

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
function _assemble!(a::Assembler{T}, dofs::AbstractVector{Int}, Ke::AbstractMatrix{T}) where {T}
    _assemble!(a, dofs, dofs, Ke)
end

""" Assembles the matrix `Ke` into `a` according to the dofs specified by `rowdofs` and `coldofs`. """
function _assemble!(a::Assembler{T}, nrows::AbstractVector{Int}, ncols::AbstractVector{Int}, Ke::AbstractMatrix{T}) where {T}
    nrows = length(rowdofs)
    ncols = length(coldofs)
    @assert size(Ke) == (nrows, ncols)

    append!(a.V, Ke)
    append!(a.I, rowdofs)
    append!(a.J, coldofs)

end