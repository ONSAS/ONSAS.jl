"Module defining a Saint-Venant-Kirchhoff material model."
module SvkMaterial

using LinearAlgebra, SparseArrays, Reexport

using ..HyperElasticMaterials
using ..Utils

@reexport import ..LinearElasticMaterials: lame_parameters, elasticity_modulus, shear_modulus,
                                           bulk_modulus, poisson_ratio
@reexport import ..HyperElasticMaterials: cosserat_stress, strain_energy

export Svk

"""
Material with Saint-Venant-Kirchhoff properties.
The strain energy `Ψ` is: `Ψ(𝔼)` = `λ`/2 tr(`𝔼`)^2 + `G` tr(`𝔼`^2).

For context see the [Hyperelastic material](https://en.wikipedia.org/wiki/Hyperelastic_material) wikipedia article.

It is also possible to construct an `Svk` material given its elasticity and shear modulus `E`, `ν` respectively and its density `ρ`.
For context see the [Lamé parameters](https://en.wikipedia.org/wiki/Lam%C3%A9_parameters) wikipedia article.
"""
struct Svk{T<:Real} <: AbstractHyperElasticMaterial
    "First Lamé parameter."
    λ::T
    "Shear modulus or second Lamé parameter (μ)."
    G::T
    "Density (`nothing` for static cases)."
    ρ::Density
    "Material label."
    label::Label
    function Svk(λ::T, G::T, ρ::Density, label::Label=NO_LABEL) where {T<:Real}
        if ρ isa Real
            ρ > 0 || error("Density must be positive.")
        end
        @assert λ ≥ 0 "The first Lamé parameter `λ` must be positive."
        @assert G ≥ 0 "The second Lamé parameter or shear modulus `G` must be positive."
        new{T}(λ, G, ρ, Symbol(label))
    end
end

"Constructor for `Svk` material with no density."
function Svk(λ::T, G::T, label::Label=NO_LABEL) where {T<:Real}
    Svk(λ, G, nothing, label)
end

"Constructor from elasticity and shear modulus `E`, `ν` respectively and density `ρ`."
function Svk(; E::Real, ν::Real, ρ::Density=nothing, label::Label=NO_LABEL)
    λ = E * ν / ((1 + ν) * (1 - 2 * ν))
    G = E / (2 * (1 + ν))
    Svk(λ, G, ρ, label)
end

"Return the strain energy for a `Svk` material `m` and the Green-Lagrange strain tensor `𝔼`."
function strain_energy(m::Svk, 𝔼::AbstractMatrix)
    λ, G = lame_parameters(m)
    λ / 2 * tr(𝔼)^2 + G * tr(𝔼^2)
end

"Return Lamé parameters `λ` and `G` from a `Svk` material `m`."
lame_parameters(m::Svk) = m.λ, m.G

"Return the shear modulus `G` from a `Svk` material `m`."
shear_modulus(m::Svk) = m.G

"Return the Poisson's ration `ν` form a `Svk` material `m`."
function poisson_ratio(m::Svk)
    λ, G = lame_parameters(m)
    λ / (2 * (λ + G))
end

"Return the elasticity modulus `E` form a `Svk` material `m`."
function elasticity_modulus(m::Svk)
    λ, G = lame_parameters(m)
    G * (3 * λ + 2 * G) / (λ + G)
end

"Return the bulk_modulus `K` for a `Svk` material `m`."
function bulk_modulus(m::Svk)
    λ, G = lame_parameters(m)
    λ + 2 * G / 3
end

"Return the Cosserat or Second-Piola Kirchoff stress tensor `𝕊`
considering a `Svk` material `m` and the Lagrangian Green
strain tensor `𝔼`.Also this function provides `∂𝕊∂𝔼` for the iterative method."
function cosserat_stress(m::Svk, 𝔼::AbstractMatrix)
    λ, G = lame_parameters(m)
    𝕊 = λ * tr(𝔼) * eye(3) + 2 * G * 𝔼

    ∂𝕊∂𝔼 = zeros(6, 6)
    ∂𝕊∂𝔼[1:3, 1:3] = λ * ones(3, 3) + 2 * G * eye(3)
    ∂𝕊∂𝔼[4:6, 4:6] = G * eye(3)

    𝕊, ∂𝕊∂𝔼
end

end
