using ..Materials: SVK
using ..Elements: AbstractElement, AbstractNode, Node
using ..CrossSections: AbstractCrossSection, area
using ..Utils: eye, row_vector

import ..Elements: local_dof_symbol

export Truss

"""
A `Truss` represents a 2D element that transmits axial force only.
### Fields:
- `n₁`             -- stores first truss node.
- `n₂`             -- stores second truss node.
- `cross_sections` -- stores the truss cross-section properties.
- `label`          -- stores the truss label.
"""
struct Truss{dim,G<:AbstractCrossSection,T<:Real} <: AbstractElement{dim,T}
    n₁::AbstractNode{dim,T}
    n₂::AbstractNode{dim,T}
    cross_section::G
    label::Symbol
    function Truss(n₁::AbstractNode{dim,T}, n₂::AbstractNode{dim,T}, g::G, label=:no_labelled_elem) where
    {dim,G<:AbstractCrossSection,T<:Real}
        new{dim,G,T}(n₁, n₂, g, Symbol(label))
    end
end

local_dof_symbol(::Truss) = [:u]
nodes(t::Truss) = [t.n₁, t.n₂]

_strain(l_ini::Number, l_def::Number) = (l_def^2 - l_ini^2) / (l_ini * (l_ini + l_def))# green lagrange strain

function strain(e::Truss, u_e::AbstractVector)
    l_def, l_ini = lengths(e, u_e)
    ϵ = _strain(l_ini, l_def)
end

function internal_forces(m::SVK, e::Truss{dim}, u_glob::AbstractVector) where {dim}

    E = m.E
    A = area(cross_section(e))

    u_e = u_glob[local_dofs(e)]

    X_ref, X_def = _X_rows(e, u_e)
    l_ref, l_def = _lengths(X_ref, X_def, dim)
    B_dif, _ = _aux_matrices(dim)

    # normalized reference and deformed co-rotational vector
    e₁_def = B_dif * X_def / l_def
    TTcl = B_dif' * e₁_def

    ϵ_e = _strain(l_ref, l_def)
    σ_e = E * ϵ_e
    fᵢₙₜ_e = A * σ_e * TTcl

    Kₘ = E * A / l_ref * (TTcl * (TTcl'))
    K_geo = σ_e * A / l_def * (B_dif' * B_dif - TTcl * (TTcl'))
    Kᵢₙₜ_e = Kₘ + K_geo

    return fᵢₙₜ_e, Kᵢₙₜ_e, σ_e, ϵ_e
end

"Returns local dofs of a truss element."
local_dofs_symbol(::Truss) = :u

function _aux_matrices(dim::Integer)
    Bdif = hcat(-eye(dim), eye(dim))
    Ge = Bdif' * Bdif
    return Bdif, Ge
end

function _X_rows(e::Truss{dim}, u_e::AbstractVector) where {dim}
    X_ref_row = reduce(vcat, coordinates(e))
    X_def_row = X_ref_row + u_e

    return X_ref_row, X_def_row
end

function _aux_b(X_ref_row::AbstractVector, X_def_row::AbstractVector, u_loc_dofs::AbstractVector, G::AbstractMatrix, dim::Integer)

    l_ref, l_def = _lengths(X_ref_row, X_def_row, dim)

    b_ref = 1 / (l_ref^2) * X_ref_row' * G
    b_def = 1 / (l_def^2) * u_loc_dofs' * G

    return b_ref, b_def
end

"Returns deformed and reference lengths of a truss element."
function _lengths(X_ref_row::AbstractVector, X_def_row::AbstractVector, dim::Integer)
    Bdif, _ = _aux_matrices(dim)
    l_ref = sqrt(sum((Bdif * X_ref_row) .^ 2))
    l_def = sqrt(sum((Bdif * X_def_row) .^ 2))

    return l_ref, l_def
end