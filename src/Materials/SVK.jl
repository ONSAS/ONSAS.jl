using LinearAlgebra: tr
using SparseArrays: SparseMatrixCSC

using .Materials: AbstractMaterial
using ..Utils: label, eye

import .Materials: density, parameters, cosserat, strain_energy, elasticity_modulus,
    shear_modulus, bulk_modulus, poisson_ratio

export SVK, lame_parameters

""" SVK material struct.
### Fields:
- `λ`     -- first Lamé parameter.
- `G`     -- shear modulus or second Lamé parameter (μ).
- `ρ`     -- density (`nothing` for static cases).
- `label` -- material label.

[See this ref.](https://en.wikipedia.org/wiki/Hyperelastic_material)
"""
struct SVK{T<:Real,R<:Union{T,Nothing}} <: AbstractMaterial
    λ::T
    G::T
    ρ::R
    label::Symbol
    function SVK(λ::T, G::T, ρ::R, label::L=:no_labelled_mat) where
    {T<:Real,R<:Union{Nothing,Real},L<:Union{Symbol,String}}
        return new{T,R}(λ, G, ρ, Symbol(label))
    end
end

"Material `SVK` constructor with no density parameter `ρ`."
function SVK(λ::Real, G::Real, label::L=:no_labelled_mat) where {L<:Union{Symbol,String}}
    return SVK(λ, G, nothing, label)
end

"Material `SVK` constructor with elasticity and shear modulus `E`, `ν` and density `ρ`. 
See [this ref](https://en.wikipedia.org/wiki/Lam%C3%A9_parameters)."
function SVK(; E::Real, ν::Real, ρ::R=nothing, label::L=:no_labelled_mat) where
{R<:Union{Nothing,Real},L<:Union{Symbol,String}}

    # Compute λ and μ (μ = G) given E and ν
    λ = E * ν / ((1 + ν) * (1 - 2 * ν))
    G = E / (2 * (1 + ν))

    return SVK(λ, G, ρ, Symbol(label))
end

"Returns the strain energy expression for a `SVK` material `m`."
strain_energy(::SVK) = :(λ / 2 * tr(𝔼)^2 + G * tr(𝔼^2))

"Returns lamé parameters `λ` and `G` from a `SVK` material `m`."
lame_parameters(m::SVK) = m.λ, m.G

"Returns the shear modulus `G` from a `SVK` material `m`."
shear_modulus(m::SVK) = m.G

"Returns the density `ρ` from a `SVK` material `m`."
density(m::SVK) = m.ρ

"Returns the Poisson's ration `ν` form a `SVK` material `m`."
function poisson_ratio(m::SVK)
    λ, G = lame_parameters(m)
    return λ / (2 * (λ + G))
end

"Returns the elasticity modulus `E` form a `SVK` material `m`."
function elasticity_modulus(m::SVK)
    λ, G = lame_parameters(m)
    return G * (3 * λ + 2 * G) / (λ + G)
end

"Returns the bulk_modulus `K` for a `SVK` material `m`."
function bulk_modulus(m::SVK)
    λ, G = lame_parameters(m)
    return λ + 2 * G / 3
end

"Returns the Cosserat or Second-Piola Kirchoff tensor (𝕊) for a `Tetrahedron` element `t`
considering a `SVK` material `m` and the Lagrangian Green strain tensor `𝔼`."
function cosserat(m::SVK, 𝔼::AbstractMatrix, compute∂𝕊∂𝔼::Bool=true)

    λ, G = lame_parameters(m)
    𝕊 = λ * tr(𝔼) * eye(3) + 2 * G * 𝔼

    if compute∂𝕊∂𝔼
        ∂𝕊∂𝔼 = SparseMatrixCSC(zeros(6, 6))
        ∂𝕊∂𝔼[1:3, 1:3] = λ * ones(3, 3) + 2 * G * eye(3)
        ∂𝕊∂𝔼[4:6, 4:6] = G * eye(3)
        return 𝕊, ∂𝕊∂𝔼
    else
        return 𝕊
    end

end


