"Module defining the interface to handle linear elastic materials."
module LinearElasticMaterials

using ..Materials

export AbstractLinearElasticMaterial, lame_parameters, elasticity_modulus, shear_modulus,
       bulk_modulus, poisson_ratio, cauchy_stress

""" Abstract supertype for all elastic material models.

An `AbstractLinearElasticMaterial` object facilitates the process of using elastic materials.

**Common methods:**
* [`elastic_modulus`](@ref)
* [`poisson_ratio`](@ref)
* [`shear_modulus`](@ref)
* [`lamé_parameters`](@ref)

**Common fields:**
* label
* ρ(density)
"""
abstract type AbstractLinearElasticMaterial <: AbstractMaterial end

"Return Lamé parameters `λ` and `G` from an `AbstractLinearElasticMaterial` material `m`."
function lame_parameters(m::AbstractLinearElasticMaterial) end

"Return shear modulus G` from an `AbstractLinearElasticMaterial` material `m`."
function shear_modulus(m::AbstractLinearElasticMaterial) end

"Return the Poisson's ratio `ν` form an `AbstractLinearElasticMaterial` material `m`."
function poisson_ratio(m::AbstractLinearElasticMaterial) end

"Return the Elasticity modulus `E` form an `AbstractLinearElasticMaterial` material `m`."
function elasticity_modulus(m::AbstractLinearElasticMaterial) end

"Return the bulk modulus `K` for an `AbstractLinearElasticMaterial` material `m`."
function bulk_modulus(m::AbstractLinearElasticMaterial) end

"Return the cauchy stress tensor `σ` and the constitutive driver `∂σ∂ϵ` 
considering a `IsotropicLinearElastic` material `m`."
function cauchy_stress(m::AbstractLinearElasticMaterial, ϵ::AbstractMatrix) end

end
