"""
Module defining the materials implemented.
Each material consists of a data type with a label and one or more parameters into its fields.
"""
module Materials

using Reexport: @reexport

@reexport import ..Utils: label

export AbstractMaterial, density, parameters
export AbstractHyperElasticMaterial, cosserat_stress, strain_energy
export AbstractLinearElasticMaterial, lame_parameters, elasticity_modulus, shear_modulus,
       bulk_modulus, poisson_ratio

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

"Return the parameters of type `Number` in the `AbstractMaterial` `m`."
function parameters(m::T) where {T<:AbstractMaterial}
    return Tuple([getfield(f, n) for n in fieldlabels(T) if fieldtype(T, n) isa Number])
end

"Return the `AbstractMaterial` `m` density `ρ`."
density(m::AbstractMaterial) = m.ρ

"Return the `AbstractMaterial` `m` label."
label(m::AbstractMaterial) = m.label

include("./LinearElasticMaterials.jl")
@reexport using .LinearElasticMaterials

include("./HyperElasticMaterials.jl")
@reexport using .HyperElasticMaterials

end # module
