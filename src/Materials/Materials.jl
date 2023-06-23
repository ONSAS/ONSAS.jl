"""
Module defining the materials implemented.
Overall, each materials consists of a data type with a label, parameters and a density into its fields.
"""
module Materials

using Reexport

@reexport import ..Utils: label

export AbstractMaterial, density, parameters

""" Abstract supertype for all material models.

An `AbstractMaterial` object facilitates the process of defining new material models.
Different material models leads to different constitutive laws, internal forces and stiffness matrices.

**Abstract Methods**

* [`parameters`](@ref)
* [`density`](@ref)
* [`label`](@ref)

**Abstract fields**
* label
* ρ(density)
"""
abstract type AbstractMaterial end

"Return the parameters of type `Number` in the `AbstractMaterial` `m`."
function parameters(m::T) where {T<:AbstractMaterial}
    Tuple([getfield(f, n) for n in fieldlabels(T) if fieldtype(T, n) isa Number])
end

"Return the `AbstractMaterial` `m` density `ρ`."
density(m::AbstractMaterial) = m.ρ

"Return the `AbstractMaterial` `m` label."
label(m::AbstractMaterial) = m.label

end
