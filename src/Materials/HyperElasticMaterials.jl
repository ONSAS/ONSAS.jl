"Module defining the interface to handle hyper-elastic materials."
module HyperElasticMaterials

using ..Materials: AbstractMaterial

export AbstractHyperElasticMaterial, cosserat_stress!, strain_energy

""" Abstract supertype for all hyper-elastic material models.

An `AbstractHyperElasticMaterial` object facilitates the process of using hyper elastic materials.
These materials are characterized by a strain energy function ψ  that depends only on the deformation gradient tensor `∇u`.

**Abstract Methods**
* [`strain_energy`](@ref)
* [`cosserat_stress!`](@ref)

**Abstract fields**
* label
* ρ(density)
"""
abstract type AbstractHyperElasticMaterial <: AbstractMaterial end

"Return the strain energy value for an `AbstractMaterial` `m`, and the Green-Lagrange
strain tensor `𝔼`."
function strain_energy(m::AbstractHyperElasticMaterial, 𝔼) end

"Return the Cosserat or Second-Piola Kirchhoff stress tensor `𝕊` given an `AbstractMaterial` `m` and the
Green-Lagrange strain tensor `𝔼`."
function cosserat_stress!(m::AbstractMaterial, 𝔼::AbstractMatrix) end

end
