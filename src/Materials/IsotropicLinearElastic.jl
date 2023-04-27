using LinearAlgebra: Symmetric, tr
using Reexport

using ..LinearElasticMaterials: AbstractLinearElasticMaterial
using ...Utils

@reexport import ..LinearElasticMaterials: lame_parameters, shear_modulus, poisson_ratio,
                                           elasticity_modulus,
                                           bulk_modulus, cauchy_stress

export IsotropicLinearElastic

"""
Material with linear elastic properties.

For context see the wikipedia article on [Linear elasticity](https://en.wikipedia.org/wiki/Linear_elasticity).

It is also possible to construct an `IsotropicLinearElastic` material given its Lamé parameters `λ`, `G` and density `ρ`. 
For context see the wikipedia article on [Lamé parameters](https://en.wikipedia.org/wiki/Lam%C3%A9_parameters).
"""
struct IsotropicLinearElastic{ET<:Real,NT<:Real} <: AbstractLinearElasticMaterial
    "Elasticity modulus."
    E::ET
    "Poisson's ratio."
    ν::NT
    "Density (`nothing` for static cases)."
    ρ::Density
    "Material label."
    label::Label
    function IsotropicLinearElastic(E::ET, ν::NT, ρ::Density,
                                    label::Label=NO_LABEL) where {ET<:Real,NT<:Real}
        new{ET,NT}(E, ν, ρ, Symbol(label))
    end
end
function IsotropicLinearElastic(E::ET, ν::NT, label::Label=NO_LABEL) where {ET<:Real,NT<:Real}
    IsotropicLinearElastic(E, ν, nothing, label)
end
function IsotropicLinearElastic(; λ::Real, G::Real, ρ::Density=nothing, label::Label=NO_LABEL)
    E = G * (3λ + 2G) / (λ + G)
    ν = λ / (2 * (λ + G))
    IsotropicLinearElastic(E, ν, ρ, Symbol(label))
end

"Return the Elasticity modulus `E` from a `IsotropicLinearElastic` material `m`."
elasticity_modulus(m::IsotropicLinearElastic) = m.E

"Return the Poisson's ratio `ν` from a `IsotropicLinearElastic` material `m`."
poisson_ratio(m::IsotropicLinearElastic) = m.ν

"Return the shear modulus `G` from a `IsotropicLinearElastic` material `m`."
shear_modulus(m::IsotropicLinearElastic) = elasticity_modulus(m) / (2 * (1 + poisson_ratio(m)))

"Return the bulk modulus `K` from a `IsotropicLinearElastic` material `m`."
bulk_modulus(m::IsotropicLinearElastic) = elasticity_modulus(m) / (3 * (1 - 2 * poisson_ratio(m)))

"Return Lamé parameters `λ` and `G` from a `IsotropicLinearElastic` material `m`."
function lame_parameters(m::IsotropicLinearElastic)
    E = elasticity_modulus(m)
    G = shear_modulus(m)
    ν = poisson_ratio(m)
    λ = E * ν / ((1 + ν) * (1 - 2 * ν))
    return λ, G
end

"Return the cauchy stress tensor `σ` and the constitutive driver `∂σ∂ϵ` 
considering a `IsotropicLinearElastic` material `m`."
function cauchy_stress(m::IsotropicLinearElastic, ϵ::AbstractMatrix)
    λ, G = lame_parameters(m)
    σ = Symmetric(λ * tr(ϵ) * eye(3) + 2 * G * ϵ)

    ∂σ∂ϵ = zeros(6, 6)
    ∂σ∂ϵ[1:3, 1:3] = λ * ones(3, 3) + 2 * G * eye(3)
    ∂σ∂ϵ[4:6, 4:6] = G * eye(3)

    return σ, ∂σ∂ϵ
end
