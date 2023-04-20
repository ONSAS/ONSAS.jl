using ForwardDiff: gradient!
using LinearAlgebra: Symmetric, tr, det, inv

using ..HyperElasticMaterials: AbstractHyperElasticMaterial
using ...Utils: eye

import ..LinearElasticMaterials: lame_parameters, elasticity_modulus, shear_modulus, bulk_modulus, poisson_ratio
import ..HyperElasticMaterials: cosserat_stress, strain_energy

export NeoHookean

""" Neo-Hookean material struct.

The strain energy `Ψ` is: `Ψ(𝔼)` = `G`/2 (tr(`ℂ`) -2 *log(`J`))^2 + `K`/2 (`J` - 1)^2

### Fields:
- `G`     -- shear modulus `G` or second Lamé parameter `μ`.
- `K`     -- bulk modulus.
- `ρ`     -- density (`nothing` for static cases).
- `label` -- material label.

[See this ref.](https://en.wikipedia.org/wiki/Neo-Hookean_solid)
"""
struct NeoHookean{T<:Real,R<:Union{T,Nothing}} <: AbstractHyperElasticMaterial
    K::T
    G::T
    ρ::R
    label::Symbol
    function NeoHookean(K::T, G::T, ρ::R, label::L=:no_labelled_mat) where
    {T<:Real,R<:Union{Nothing,Real},L<:Union{Symbol,String}}
        if ρ isa Real
            ρ > 0 || error("Density must be positive.")
        end
        @assert K ≥ 0 "The bulk modulus `K` must be positive."
        @assert G ≥ 0 "The shear modulus or second Lamé parameter `μ` must be positive."
        return new{T,R}(K, G, ρ, Symbol(label))
    end
end

"Material `NeoHookean` constructor with no density parameter `ρ`."
function NeoHookean(K::Real, G::Real, label::L=:no_labelled_mat) where {L<:Union{Symbol,String}}
    return NeoHookean(K, G, nothing, label)
end

"Material `NeoHookean` constructor with elasticity and shear modulus `E`, `ν` and density `ρ`. 
See [this ref](https://en.wikipedia.org/wiki/Lam%C3%A9_parameters)."
function NeoHookean(; E::Real, ν::Real, ρ::R=nothing, label::L=:no_labelled_mat) where
{R<:Union{Nothing,Real},L<:Union{Symbol,String}}

    # Compute λ, μ and K (μ = G) given E and ν
    λ = E * ν / ((1 + ν) * (1 - 2 * ν))
    G = E / (2 * (1 + ν))
    K = λ + 2 * G / 3

    return NeoHookean(K, G, ρ, Symbol(label))
end
"Return the strain energy for a `NeoHookean` material `m` and the Green-Lagrange strain tensor `𝔼`."
function strain_energy(m::NeoHookean, 𝔼::AbstractMatrix)
    ℂ = Symmetric(2 * 𝔼 + eye(3))
    J = sqrt(det(ℂ))
    # First invariant
    I₁ = tr(ℂ)
    # Strain energy function 
    Ψ = shear_modulus(m) / 2 * (I₁ - 2 * log(J)) + bulk_modulus(m) / 2 * (J - 1)^2
end

"Return Lamé parameters `λ` and `G` from a `NeoHookean` material `m`."
function lame_parameters(m::NeoHookean)
    G = shear_modulus(m)
    λ = bulk_modulus(m) - 2 * G / 3
    return λ, G
end

"Return the shear modulus `G` from a `NeoHookean` material `m`."
shear_modulus(m::NeoHookean) = m.G

"Return the Poisson's ration `ν` form a `NeoHookean` material `m`."
function poisson_ratio(m::NeoHookean)
    λ, G = lame_parameters(m)
    λ / (2 * (λ + G))
end

"Return the elasticity modulus `E` form a `NeoHookean` material `m`."
function elasticity_modulus(m::NeoHookean)
    λ, G = lame_parameters(m)
    G * (3 * λ + 2 * G) / (λ + G)
end

"Return the bulk_modulus `K` for a `NeoHookean` material `m`."
bulk_modulus(m::NeoHookean) = m.K

"Return the Cosserat stress tensor `𝕊` given the Green-Lagrange `𝔼` strain tensor."
function _𝕊_analytic(m::NeoHookean, 𝔼::AbstractMatrix)
    # Right hand Cauchy strain tensor
    ℂ = Symmetric(2 * 𝔼 + eye(3))
    ℂ⁻¹ = inv(ℂ)
    J = sqrt(det(ℂ))
    # Compute 𝕊 
    shear_modulus(m) * (eye(3) - ℂ⁻¹) + bulk_modulus(m) * (J * (J - 1) * ℂ⁻¹)
end

"Return the `∂𝕊∂𝔼` for a material `m`, the Gree-Lagrange strain tensor `𝔼` and a
function to compute 𝕊 analytically."
function _∂𝕊_∂𝔼(m::NeoHookean, 𝔼::AbstractMatrix, 𝕊_analytic::Function=_𝕊_analytic)

    indexes = [(1, 1), (2, 2), (3, 3), (2, 3), (1, 3), (1, 2)]

    ∂S∂𝔼_forward_diff = zeros(6, 6)
    aux_gradients = zeros(3, 3)

    row = 1
    for index in indexes
        i, j = index
        ∂S∂𝔼_forward_diff[row, :] .= _voigt(
            gradient!(aux_gradients, E -> 𝕊_analytic(m, E)[i, j], collect(𝔼)), #TODO: Fix with Symmetric
            0.5
        )
        row += 1
    end
    return ∂S∂𝔼_forward_diff

end

"Return the Cosserat or Second-Piola Kirchoff stress tensor `𝕊` 
considering a `SVK` material `m` and the Green-Lagrange  
strain tensor `𝔼`.Also this function provides `∂𝕊∂𝔼` for the iterative method."
function cosserat_stress(m::NeoHookean, 𝔼::AbstractMatrix)

    return _𝕊_analytic(m, 𝔼), _∂𝕊_∂𝔼(m, 𝔼, _𝕊_analytic)

end
