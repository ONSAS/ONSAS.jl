using LinearAlgebra: tr
using SparseArrays: SparseMatrixCSC
using Reexport

using ..HyperElasticMaterials
using ...Utils

@reexport import ..LinearElasticMaterials: lame_parameters, elasticity_modulus, shear_modulus,
                                           bulk_modulus,
                                           poisson_ratio
@reexport import ..HyperElasticMaterials: cosserat_stress, strain_energy

export SVK

""" 
Material with Saint-Venant-Kirchhoff properties.
The strain energy `Ψ` is: `Ψ(𝔼)` = `λ`/2 tr(`𝔼`)^2 + `G` tr(`𝔼`^2).

For context see the [Hyperelastic material](https://en.wikipedia.org/wiki/Hyperelastic_material) wikipedia article.

It is also possible to construct an `SVK` material given its elasticity and shear modulus `E`, `ν` respectively and its density `ρ`. 
For context see the [Lamé parameters](https://en.wikipedia.org/wiki/Lam%C3%A9_parameters) wikipedia article.
"""
struct SVK{T<:Real,R<:Union{T,Nothing}} <: AbstractHyperElasticMaterial
    "First Lamé parameter."
    λ::T
    "Shear modulus or second Lamé parameter (μ)."
    G::T
    "Density (`nothing` for static cases)."
    ρ::R
    "Material label."
    label::Label
    function SVK(λ::T, G::T, ρ::R, label::Label=NO_LABEL) where {T<:Real,R<:Union{Nothing,Real}}
        if ρ isa Real
            ρ > 0 || error("Density must be positive.")
        end
        @assert λ ≥ 0 "The first Lamé parameter `λ` must be positive."
        @assert G ≥ 0 "The second Lamé parameter or shear modulus `G` must be positive."
        return new{T,R}(λ, G, ρ, Symbol(label))
    end
end
function SVK(λ::T, G::T, label::Label=NO_LABEL) where {T<:Real,R<:Union{Nothing,Real}}
    SVK(λ, G, nothing, label)
end
function SVK(; E::Real, ν::Real, ρ::R=nothing, label::Label=NO_LABEL) where {R<:Union{Nothing,Real}}
    # Compute λ and μ (μ = G) given E and ν.
    λ = E * ν / ((1 + ν) * (1 - 2 * ν))
    G = E / (2 * (1 + ν))
    return SVK(λ, G, ρ, label)
end

"Return the strain energy for a `SVK` material `m` and the Green-Lagrange strain tensor `𝔼`."
function strain_energy(m::SVK, 𝔼::AbstractMatrix)
    λ, G = lame_parameters(m)
    return λ / 2 * tr(𝔼)^2 + G * tr(𝔼^2)
end

"Return Lamé parameters `λ` and `G` from a `SVK` material `m`."
lame_parameters(m::SVK) = m.λ, m.G

"Return the shear modulus `G` from a `SVK` material `m`."
shear_modulus(m::SVK) = m.G

"Return the Poisson's ration `ν` form a `SVK` material `m`."
function poisson_ratio(m::SVK)
    λ, G = lame_parameters(m)
    return λ / (2 * (λ + G))
end

"Return the elasticity modulus `E` form a `SVK` material `m`."
function elasticity_modulus(m::SVK)
    λ, G = lame_parameters(m)
    return G * (3 * λ + 2 * G) / (λ + G)
end

"Return the bulk_modulus `K` for a `SVK` material `m`."
function bulk_modulus(m::SVK)
    λ, G = lame_parameters(m)
    return λ + 2 * G / 3
end

"Return the Cosserat or Second-Piola Kirchoff stress tensor `𝕊` 
considering a `SVK` material `m` and the Lagrangian Green 
strain tensor `𝔼`.Also this function provides `∂𝕊∂𝔼` for the iterative method."
function cosserat_stress(m::SVK, 𝔼::AbstractMatrix)
    λ, G = lame_parameters(m)
    𝕊 = λ * tr(𝔼) * eye(3) + 2 * G * 𝔼

    ∂𝕊∂𝔼 = zeros(6, 6)
    ∂𝕊∂𝔼[1:3, 1:3] = λ * ones(3, 3) + 2 * G * eye(3)
    ∂𝕊∂𝔼[4:6, 4:6] = G * eye(3)

    return 𝕊, ∂𝕊∂𝔼
end
