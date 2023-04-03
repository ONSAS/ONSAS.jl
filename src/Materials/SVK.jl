using LinearAlgebra: tr
using SparseArrays: SparseMatrixCSC

using .Materials: AbstractMaterial
using ..Utils: eye

import .Materials: density, cosserat, strain_energy, lame_parameters, elasticity_modulus, shear_modulus, bulk_modulus, poisson_ratio
export SVK

""" Saint-Venant-Kirchhoff material struct.

The strain energy `Ψ` is: `Ψ(𝔼)` = `λ`/2 tr(`𝔼`)^2 + `G` tr(`𝔼`^2)

### Fields:
- `λ`     -- first Lamé parameter.
- `G`     -- shear modulus or second Lamé parameter (μ).
- `ρ`     -- density (`nothing` for static cases).
- `label` -- material label.

[See this ref.](https://en.wikipedia.org/wiki/Hyperelastic_material)
"""
struct SVK{T<:Real,R<:Union{T,Nothing}} <: AbstractHyperElasticMaterial
    λ::T
    G::T
    ρ::R
    label::Symbol
    function SVK(λ::T, G::T, ρ::R, label::L=:no_labelled_mat) where
    {T<:Real,R<:Union{Nothing,Real},L<:Union{Symbol,String}}
        if ρ isa Real
            ρ > 0 || error("Density must be positive.")
        end
        @assert λ ≥ 0 "The first Lamé parameter `λ` must be positive."
        @assert G ≥ 0 "The second Lamé parameter or shear modulus `G` must be positive."
        return new{T,R}(λ, G, ρ, Symbol(label))
    end
end


"Material `SVK` constructor with no density parameter `ρ`."
SVK(λ::Real, G::Real, label::L=:no_labelled_mat) where {L<:Union{Symbol,String}} =
    SVK(λ, G, nothing, label)

"Material `SVK` constructor with elasticity and shear modulus `E`, `ν` and density `ρ`. 
See [this ref](https://en.wikipedia.org/wiki/Lam%C3%A9_parameters)."
function SVK(; E::Real, ν::Real, ρ::R=nothing, label::L=:no_labelled_mat) where
{R<:Union{Nothing,Real},L<:Union{Symbol,String}}

    # Compute λ and μ (μ = G) given E and ν
    λ = E * ν / ((1 + ν) * (1 - 2 * ν))
    G = E / (2 * (1 + ν))

    SVK(λ, G, ρ, Symbol(label))
end

"Returns the strain energy for a `SVK` material `m` and the Green-Lagrange strain tensor `𝔼`."
function strain_energy(m::SVK, 𝔼::AbstractMatrix)
    λ, G = lame_parameters(m)
    λ / 2 * tr(𝔼)^2 + G * tr(𝔼^2)
end

"Returns Lamé parameters `λ` and `G` from a `SVK` material `m`."
lame_parameters(m::SVK) = m.λ, m.G

"Returns the shear modulus `G` from a `SVK` material `m`."
shear_modulus(m::SVK) = m.G

"Returns the Poisson's ration `ν` form a `SVK` material `m`."
function poisson_ratio(m::SVK)
    λ, G = lame_parameters(m)
    λ / (2 * (λ + G))
end

"Returns the elasticity modulus `E` form a `SVK` material `m`."
function elasticity_modulus(m::SVK)
    λ, G = lame_parameters(m)
    G * (3 * λ + 2 * G) / (λ + G)
end

"Returns the bulk_modulus `K` for a `SVK` material `m`."
function bulk_modulus(m::SVK)
    λ, G = lame_parameters(m)
    λ + 2 * G / 3
end

"Returns the Cosserat or Second-Piola Kirchoff stress tensor `𝕊` 
considering a `SVK` material `m` and the Lagrangian Green 
strain tensor `𝔼`.Also this function provides `∂𝕊∂𝔼` for the iterative method."
function cosserat(m::SVK, 𝔼::AbstractMatrix)

    λ, G = lame_parameters(m)
    𝕊 = λ * tr(𝔼) * eye(3) + 2 * G * 𝔼

    ∂𝕊∂𝔼 = SparseMatrixCSC(zeros(6, 6))
    ∂𝕊∂𝔼[1:3, 1:3] = λ * ones(3, 3) + 2 * G * eye(3)
    ∂𝕊∂𝔼[4:6, 4:6] = G * eye(3)

    return 𝕊, ∂𝕊∂𝔼

end