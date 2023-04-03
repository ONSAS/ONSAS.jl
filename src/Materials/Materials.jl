"""
Module defining the materials implemented.
Each material consists of a data type with a label and one or more parameters into its fields.
"""
module Materials

using Reexport: @reexport

@reexport import ..Utils: label

export AbstractMaterial, density, parameters
export AbstractHyperElasticMaterial, cosserat, strain_energy
export AbstractLinearElasticMaterial, lame_parameters, elasticity_modulus, shear_modulus, bulk_modulus, poisson_ratio

""" Abstract supertype for all material models.

An `AbstractMaterial` object facilitates the process of defining new material models. 
Different material models leads to different constitutive laws, internal forces and stiffness matrices.

**Common methods:**

* [`parameters`](@ref)
* [`density`](@ref)
* [`label`](@ref)

**Common fields:**
* label
* ρ(density)
"""
abstract type AbstractMaterial end

"Returns the parameters of type `Number` in the `AbstractMaterial` `m`."
parameters(m::T) where {T<:AbstractMaterial} = Tuple([getfield(f, n) for n in fieldlabels(T) if fieldtype(T, n) isa Number])

"Returns the `AbstractMaterial` `m` density `ρ`."
density(m::AbstractMaterial) = m.ρ

"Returns the `AbstractMaterial` `m` label."
label(m::AbstractMaterial) = m.label

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

"Returns Lamé parameters `λ` and `G` from an `AbstractLinearElasticMaterial` material `m`."
function lame_parameters(m::AbstractLinearElasticMaterial) end

"Returns shear modulus G` from an `AbstractLinearElasticMaterial` material `m`."
function shear_modulus(m::AbstractLinearElasticMaterial) end

"Returns the Poisson's ratio `ν` form an `AbstractLinearElasticMaterial` material `m`."
function poisson_ratio(m::AbstractLinearElasticMaterial) end

"Returns the Elasticity modulus `E` form an `AbstractLinearElasticMaterial` material `m`."
function elasticity_modulus(m::AbstractLinearElasticMaterial) end

"Returns the bulk modulus `K` for an `AbstractLinearElasticMaterial` material `m`."
function bulk_modulus(m::AbstractLinearElasticMaterial) end

include("IsotropicLinearElastic.jl")

""" Abstract supertype for all hyper-elastic material models.

An `AbstractHyperElasticMaterial` object facilitates the process of using hyper elastic materials.
These materials are characterized by a strain energy function ψ  that depends only on the deformation gradient tensor `∇u`.

**Common methods:**
* [`strain_energy`](@ref)
* [`cosserat`](@ref)

**Common fields:**
* label
* ρ(density)

"""
abstract type AbstractHyperElasticMaterial <: AbstractMaterial end

"Returns the strain energy value for an `AbstractMaterial` `m`, and the Green-Lagrange 
strain tensor `𝔼`."
function strain_energy(m::AbstractHyperElasticMaterial, 𝔼) end

"Returns the Cosserat or Second-Piola Kirchhoff stress tensor `𝕊` given an `AbstractMaterial` `m` and the 
Green-Lagrange strain tensor `𝔼`."
function cosserat(m::AbstractMaterial, 𝔼::AbstractMatrix) end

include("SVK.jl")
include("NeoHookean.jl")
include("HyperElastic.jl")


end # module