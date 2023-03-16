"""
Module defining the materials implemented.
Each material consists of a data type with a label and one or more parameters into its fields.
"""
module Materials

using Reexport: @reexport

@reexport import ..Utils: label

export AbstractMaterial, density, parameters, cosserat, strain_energy

""" Abstract supertype for all material models.

An `AbstractMaterial` object facilitates the process of defining new material models. 
Different material models leads to different constitutive laws, internal forces and stiffness matrices.

**Common methods:**

* [`parameters`](@ref)
* [`density`](@ref)
* [`label`](@ref)

**Common methods for hyper-elastic materials:**
* [`cosserat`](@ref)
* [`strain_energy`](@ref)

**Common fields:**
* label
* ρ(density)
"""
abstract type AbstractMaterial end

"Returns the parameters of type `Number` in the `AbstractMaterial` `m`."
parameters(m::T) where {T<:AbstractMaterial} = Tuple([getfield(f, n) for n in fieldlabels(T) if fieldtype(T, n) isa Number])

"Returns the Cosserat or Second-Piola Kirchhoff stress tensor `𝕊` given an `AbstractMaterial` `m` and the 
Green-Lagrange strain tensor `𝔼`."
function cosserat(m::AbstractMaterial, 𝔼::AbstractMatrix) end

"Returns the `AbstractMaterial` `m` density `ρ`."
density(m::AbstractMaterial) = m.ρ

"Returns the `AbstractMaterial` `m` label."
label(m::AbstractMaterial) = m.label

"Returns the strain energy value for an `AbstractMaterial` `m`, and the Green-Lagrange 
strain tensor `𝔼`."
function strain_energy(m::AbstractMaterial, 𝔼) end

#===================================#
# AbstractMaterial implementations #
#==================================#

include("./SVK.jl")
include("./NeoHookean.jl")
include("./HyperElastic.jl")


end # module