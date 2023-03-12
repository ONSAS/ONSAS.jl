"""
Module defining the materials implemented.
Each material consists of a data type with a label and one or more parameters into its fields.
"""
module Materials

using Reexport: @reexport

@reexport import ..Utils: label

export AbstractMaterial, density, parameters, cosserat, strain_energy,
    elasticity_modulus, bulk_modulus, shear_modulus, poisson_ratio

""" Abstract supertype for all material models.

An `AbstractMaterial` object facilitates the process of defining new material models. 
Different material models leads to different constitutive laws, internal forces and stiffness matrices.

**Common methods:**


* [`strain_energy`](@ref)
* [`parameters`](@ref)
* [`density`](@ref)
* [`bulk_modulus`](@ref)
* [`elasticity_modulus`](@ref)
* [`shear_modulus`](@ref)
* [`poisson_ratio`](@ref)
* [`cosserat`](@ref)
* [`label`](@ref)
"""
abstract type AbstractMaterial end

"Returns the parameters of type `Number` in the `AbstractMaterial` `m`."
parameters(m::T) where {T<:AbstractMaterial} = Tuple([getfield(f, n) for n in fieldlabels(T) if fieldtype(T, n) isa Number])

"Returns the Cosserat or Second-Piola Kirchhoff stress tensor `𝕊``."
function cosserat(m::AbstractMaterial) end

"Returns the `AbstractMaterial` `m` density `ρ`."
function density(m::AbstractMaterial) end

"Returns the equivalent elasticity modulus `E` for the `AbstractMaterial` `m`."
function elasticity_modulus(m::AbstractMaterial) end

"Returns the equivalent Shear modulus `G` for the `AbstractMaterial` `m`."
function shear_modulus(m::AbstractMaterial) end

"Returns the equivalent Poisson's ratio `ν` for the `AbstractMaterial` `m`."
function poisson_ratio(m::AbstractMaterial) end

"Returns the equivalent Bulk modulus `K` for the `AbstractMaterial` `m`."
function bulk_modulus(m::AbstractMaterial) end

"Returns the `AbstractMaterial` `m` label."
label(m::AbstractMaterial) = m.label

"Returns the `AbstractMaterial` `m` strain energy expression ϕ(𝔼) ."
function strain_energy(m::AbstractMaterial) end

#===================================#
# AbstractMaterial implementations #
#==================================#

include("./SVK.jl")
# include("./NeoHookean.jl")

end # module