using SparseArrays: sparse
using ..Materials
using ..IsotropicLinearElasticMaterial
using ..Elements
using ..CrossSections
using ..Utils

import ..Elements: nodes, create_entity, cross_section, internal_forces, local_dof_symbol, strain,
                   stress

export Truss

"""
A `Truss` represents an element composed by two `Node`s that transmits axial force only.

### References
See [[ANLE]](@ref).
"""
struct Truss{dim,T<:Real,N<:AbstractNode{dim,T},VN<:AbstractVector{N},G<:AbstractCrossSection} <:
       AbstractElement{dim,T}
    "Stores the truss nodes."
    nodes::VN
    "Stores the truss cross-section properties."
    cross_section::G
    "Stores the truss label."
    label::Label
    function Truss(nodes::VN, g::G,
                   label::Label=NO_LABEL) where
             {dim,T<:Real,N<:AbstractNode{dim,T},VN<:AbstractVector{N},G<:AbstractCrossSection}
        @assert 1 ≤ dim ≤ 3 "Nodes of a truss element must comply  1 < dim < 3 ."
        new{dim,T,N,VN,G}(nodes, g, Symbol(label))
    end
end

"Constructor for a `Truss` element considering the nodes `n₁` and `n₂` and the cross-section `g`."
function Truss(n₁::N, n₂::N, g::G,
               label::Label=NO_LABEL) where
         {dim,T<:Real,N<:AbstractNode{dim,T},G<:AbstractCrossSection}
    Truss(SVector(n₁, n₂), g, label)
end

"Constructor for a `Truss` element without nodes and a `label`. This function is used to create meshes via GMSH."
function Truss(g::AbstractCrossSection, label::Label=NO_LABEL)
    return Truss(Node(0, 0, 0), Node(0, 0, 0), g, label)
end

#==============================#
# Truss element hard contracts #
#==============================#

"Return the cross-section of a `Truss` element `t`."
cross_section(t::Truss) = t.cross_section

"Return a `Tetrahedron` given an empty `Tetrahedron` `t` and a `Vector` of `Node`s `vn`."
function create_entity(t::Truss, vn::AbstractVector{<:AbstractNode})
    return Truss(vn[1], vn[2], cross_section(t), label(t))
end

"Return the local dof symbol of a `Truss` element."
local_dof_symbol(::Truss) = [:u]

"Return the internal force of a `Truss` element `t` formed by an `AbstractMaterial` `m` 
and a an element displacement vector `u_e`."
function internal_forces(m::AbstractMaterial, e::Truss{dim}, u_e::AbstractVector) where {dim}
    E = elasticity_modulus(m)
    A = area(cross_section(e))
    X_ref, X_def = _X_rows(e, u_e)
    l_ref, l_def = _lengths(X_ref, X_def, dim)
    B_dif, _ = _aux_matrices(dim)

    # normalized reference and deformed co-rotational vector
    e₁_def = B_dif * X_def / l_def
    TTcl = B_dif' * e₁_def

    ϵ = _strain(l_ref, l_def)
    σ = E * ϵ
    fᵢₙₜ_e = A * σ * TTcl

    Kₘ = E * A / l_ref * (TTcl * (TTcl'))
    K_geo = σ * A / l_def * (B_dif' * B_dif - TTcl * (TTcl'))
    Kᵢₙₜ_e = Kₘ + K_geo

    σ_e = sparse(zeros(3, 3))
    ϵ_e = sparse(zeros(3, 3))
    σ_e[1, 1] = σ
    ϵ_e[1, 1] = ϵ

    return fᵢₙₜ_e, Kᵢₙₜ_e, σ_e, ϵ_e
end

"Return the rotated engineering strain for a given reference and deformed length `l_ini` and the deformed length `l_def`. "
_strain(l_ini::Number, l_def::Number) = (l_def^2 - l_ini^2) / (l_ini * (l_ini + l_def))# rotated engi lagrange strain

"Return the strain of given `Truss` element `t` with a element displacement vector `u_e`. "
function strain(t::Truss{dim}, u_e::AbstractVector) where {dim}
    X_ref, X_def = _X_rows(t, u_e)
    l_ref, l_def = _lengths(X_ref, X_def, dim)
    return ϵ = _strain(l_ref, l_def)
end

"Return the stress of given `Truss` element `t` with a element displacement vector `u_e`. "
stress(m::IsotropicLinearElastic, t::Truss, u_e::AbstractVector) = m.E * strain(t, u_e)

"Return "
function _aux_matrices(dim::Integer)
    Bdif = hcat(-eye(dim), eye(dim))
    Ge = Bdif' * Bdif
    return Bdif, Ge
end

"Return auxiliar matrices with the element coordinates at the deformed and reference configurations."
function _X_rows(e::Truss{dim}, u_e::AbstractVector) where {dim}
    X_ref_row = reduce(vcat, coordinates(e))
    X_def_row = X_ref_row + u_e

    return X_ref_row, X_def_row
end

"Return auxiliar vectors b_ref and b_def of a truss element."
function _aux_b(X_ref_row::AbstractVector, X_def_row::AbstractVector, u_loc_dofs::AbstractVector,
                G::AbstractMatrix, dim::Integer)
    l_ref, l_def = _lengths(X_ref_row, X_def_row, dim)

    b_ref = 1 / (l_ref^2) * X_ref_row' * G
    b_def = 1 / (l_def^2) * u_loc_dofs' * G

    return b_ref, b_def
end

"Return deformed and reference lengths of a truss element."
function _lengths(X_ref_row::AbstractVector, X_def_row::AbstractVector, dim::Integer)
    Bdif, _ = _aux_matrices(dim)
    l_ref = sqrt(sum((Bdif * X_ref_row) .^ 2))
    l_def = sqrt(sum((Bdif * X_def_row) .^ 2))

    return l_ref, l_def
end
