"Module defining truss elements."
module Trusses

using SparseArrays
using Reexport
using StaticArrays

using ..Materials
using ..IsotropicLinearElasticMaterial
using ..Nodes
using ..Entities
using ..CrossSections
using ..Utils

@reexport import ..Entities: nodes, create_entity, cross_section, internal_forces, local_dof_symbol,
                             strain, stress

export Truss, strain_model
export AbstractStrainModel, RotatedEngineeringStrain, GreenStrain

"""
A `AbstractStrainModel` represents the strain model used.

### References
See [[ANLE]](@ref).
"""
abstract type AbstractStrainModel end

"""
A `RotatedEngineeringStrain` a singleton to use  the engineering strain model.

The strain is: `ϵᵢ = (lₙ - l₀) / l₀`
"""
struct RotatedEngineeringStrain <: AbstractStrainModel end

"""
A `GreenStrain` a singleton to use  the Green engineering strain model.

The strain is: `ϵᵢ = (lₙ² - l₀²) / (2l₀²)`
"""
struct GreenStrain <: AbstractStrainModel end

const DEFAULT_STRAIN_MODEL = RotatedEngineeringStrain

"""
A `Truss` represents an element composed by two `Node`s that transmits axial force only.

### References
See [[ANLE]](@ref).
"""
struct Truss{dim,E<:AbstractStrainModel,T<:Real,N<:AbstractNode{dim,T},VN<:AbstractVector{N},
             G<:AbstractCrossSection} <:
       AbstractElement{dim,T}
    "Stores the truss nodes."
    nodes::VN
    "Stores the truss cross-section properties."
    cross_section::G
    "Stores the truss label."
    label::Label
    function Truss(nodes::VN, g::G,
                   ::Type{E}=DEFAULT_STRAIN_MODEL,
                   label::Label=NO_LABEL) where
             {dim,E<:AbstractStrainModel,T<:Real,N<:AbstractNode{dim,T},VN<:AbstractVector{N},
              G<:AbstractCrossSection}
        @assert 1 ≤ dim ≤ 3 "Nodes of a truss element must comply  1 < dim < 3 ."
        new{dim,E,T,N,VN,G}(nodes, g, Symbol(label))
    end
end

"Constructor for a `Truss` element considering the nodes `n₁` and `n₂` and the cross-section `g` and 
strain model."
function Truss(n₁::N, n₂::N, g::G,
               strain::Type{E}=DEFAULT_STRAIN_MODEL,
               label::Label=NO_LABEL) where
         {dim,E<:AbstractStrainModel,T<:Real,N<:AbstractNode{dim,T},G<:AbstractCrossSection}
    Truss(SVector(n₁, n₂), g, strain, label)
end

"Constructor for a `Truss` element considering the nodes `n₁` and `n₂` and the cross-section `g`."
function Truss(n₁::N, n₂::N, g::G,
               label::Label) where
         {dim,T<:Real,N<:AbstractNode{dim,T},G<:AbstractCrossSection}
    Truss(SVector(n₁, n₂), g, DEFAULT_STRAIN_MODEL, label)
end

"Constructor for a `Truss` element without nodes a `label` and `strain`. This function is used to create meshes via GMSH."
function Truss(g::AbstractCrossSection, ::Type{E}=DEFAULT_STRAIN_MODEL,
               label::Label=NO_LABEL) where {E<:AbstractStrainModel}
    Truss(Node(0, 0, 0), Node(0, 0, 0), g, E, label)
end

"Constructor for a `Truss` element without nodes. This function is used to create meshes via GMSH."
function Truss(g::AbstractCrossSection, label::Label=NO_LABEL)
    Truss(Node(0, 0, 0), Node(0, 0, 0), g, DEFAULT_STRAIN_MODEL, label)
end

#==============================#
# Truss element hard contracts #
#==============================#

"Return the cross-section of a `Truss` element `t`."
cross_section(t::Truss) = t.cross_section

"Return the strain model used"
strain_model(::Truss{dim,E}) where {dim,E<:AbstractStrainModel} = E

"Return a `Tetrahedron` given an empty `Tetrahedron` `t` and a `Vector` of `Node`s `vn`."
function create_entity(t::Truss, vn::AbstractVector{<:AbstractNode})
    Truss(vn[1], vn[2], cross_section(t), strain_model(t), label(t))
end

"Return the local dof symbol of a `Truss` element."
local_dof_symbol(::Truss) = [:u]

"Return the internal force of a `Truss` element `t` formed by an `AbstractMaterial` `m` 
and an element displacement vector `u_e`."
function internal_forces(m::AbstractMaterial, e::Truss{dim,EM},
                         u_e::AbstractVector) where {dim,EM<:AbstractStrainModel}
    E = elasticity_modulus(m)
    A = area(cross_section(e))
    X_ref, X_def = _X_rows(e, u_e)
    l_ref, l_def = _lengths(X_ref, X_def, dim)
    B_dif, _ = _aux_matrices(dim)

    # normalized reference and deformed co-rotational vector
    e₁_def = B_dif * X_def / l_def
    TTcl = B_dif' * e₁_def

    ϵ = _strain(l_ref, l_def, EM)
    σ = E * ϵ
    fᵢₙₜ_e = A * σ * TTcl

    Kₘ = E * A / l_ref * (TTcl * (TTcl'))
    K_geo = σ * A / l_def * (B_dif' * B_dif - TTcl * (TTcl'))
    Kᵢₙₜ_e = Kₘ + K_geo

    σ_e = sparse(zeros(3, 3))
    ϵ_e = sparse(zeros(3, 3))
    σ_e[1, 1] = σ
    ϵ_e[1, 1] = ϵ

    fᵢₙₜ_e, Kᵢₙₜ_e, σ_e, ϵ_e
end

"Return the `RotatedEngineeringStrain` for a given reference and deformed length `l_ini` and the deformed length `l_def`. "
function _strain(l_ini::Real, l_def::Real, ::Type{RotatedEngineeringStrain})
    (l_def^2 - l_ini^2) / (l_ini * (l_ini + l_def))
end

"Return the `GreenStrain` for a given reference and deformed length `l_ini` and the deformed length `l_def`. "
function _strain(l_ini::Real, l_def::Real, ::Type{GreenStrain})
    (l_def^2 - l_ini^2) / (2 * l_ini^2)
end

"Return the strain of given `Truss` element `t` with a element displacement vector `u_e`. "
function strain(t::Truss{dim,E}, u_e::AbstractVector) where {dim,E<:AbstractStrainModel}
    X_ref, X_def = _X_rows(t, u_e)
    l_ref, l_def = _lengths(X_ref, X_def, dim)
    _strain(l_ref, l_def, E)
end

"Return the stress of given `Truss` element `t` with a element displacement vector `u_e`. "
stress(m::IsotropicLinearElastic, t::Truss, u_e::AbstractVector) = m.E * strain(t, u_e)

"Return "
function _aux_matrices(dim::Integer)
    Bdif = hcat(-eye(dim), eye(dim))
    Ge = Bdif' * Bdif
    Bdif, Ge
end

"Return auxiliar matrices with the element coordinates at the deformed and reference configurations."
function _X_rows(e::Truss{dim}, u_e::AbstractVector) where {dim}
    X_ref_row = reduce(vcat, coordinates(e))
    X_def_row = X_ref_row + u_e
    X_ref_row, X_def_row
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
    l_ref, l_def
end

end
