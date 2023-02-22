using ..Materials: SVK
using ..Elements: AbstractElement, AbstractNode
using ..CrossSections: AbstractCrossSection, area
using ..Utils: eye

import ..Elements: nodes, internal_forces, local_dof_symbol, strain, stress

export Truss

"""
A `Truss` represents an element composed by two `Node`s that transmits axial force only.
### Fields:
- `n₁`             -- stores first truss node.
- `n₂`             -- stores second truss node.
- `cross_sections` -- stores the truss cross-section properties.
- `label`          -- stores the truss label.
"""
struct Truss{dim,N<:AbstractNode{dim},G<:AbstractCrossSection,T<:Real} <: AbstractElement{dim,T}
    n₁::N
    n₂::N
    cross_section::G
    label::Symbol
    function Truss(n₁::N, n₂::N, g::G, label=:no_labelled_elem) where
    {dim,T<:Real,N<:AbstractNode{dim,T},G<:AbstractCrossSection}
        new{dim,N,G,T}(n₁, n₂, g, Symbol(label))
    end
end

#==============================#
# Truss element hard contracts #
#==============================#

"Returns the local dof symbol of a `Truss` element."
local_dof_symbol(::Truss) = [:u]

"Returns a vector of `Node`'s corresponding to a `Truss` element `t`."
nodes(t::Truss) = [t.n₁, t.n₂]

"Returns the internal force of a `Truss` element `t` doted with an `SVK` material and a global displacement vector `u_glob`."
function internal_forces(m::SVK, e::Truss{dim}, u_e::AbstractVector) where {dim}

    E = m.E
    A = area(cross_section(e))
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

"Returns the Green strain of given the reference length `l_ini` and the deformed length `l_def`. "
_strain(l_ini::Number, l_def::Number) = (l_def^2 - l_ini^2) / (l_ini * (l_ini + l_def))# rotated engi lagrange strain

"Returns the strain of given `Truss` element `t` with a element displacement vector `u_e`. "
function strain(t::Truss{dim}, u_e::AbstractVector) where {dim}
    X_ref, X_def = _X_rows(t, u_e)
    l_ref, l_def = _lengths(X_ref, X_def, dim)
    ϵ = _strain(l_ref, l_def)
end

"Returns the stress of given `Truss` element `t` with a element displacement vector `u_e`. "
stress(m::SVK, t::Truss, u_e::AbstractVector) = m.E * strain(t, u_e)

"Returns "
function _aux_matrices(dim::Integer)
    Bdif = hcat(-eye(dim), eye(dim))
    Ge = Bdif' * Bdif
    return Bdif, Ge
end

"Returns auxiliar matrices with the element coordinates at the deformed and reference configurations."
function _X_rows(e::Truss{dim}, u_e::AbstractVector) where {dim}
    X_ref_row = reduce(vcat, coordinates(e))
    X_def_row = X_ref_row + u_e

    return X_ref_row, X_def_row
end

"Returns auxiliar vectors b_ref and b_def of a truss element."
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