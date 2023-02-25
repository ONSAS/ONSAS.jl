"""
Module defining the materials implemented.
Each material consists of a data type with a label and one or more parameters into its fields.
"""
module Materials

using Reexport: @reexport

@reexport import ..Utils: label

export AbstractMaterial, parameters, cosserat

""" Abstract supertype for all material models.

An `AbstractMaterial` object facilitates the process of defining new material models. 
Different material models leads to different constitutive laws, internal forces and stiffness matrices.

**Common methods:**

* [`cosserat`](@ref)
* [`parameters`](@ref)
* [`label`](@ref)
"""
abstract type AbstractMaterial end

"Returns the Cosserat or Second-Piola Kirchhoff material parameters."
function cosserat(m::AbstractMaterial) end

"Returns the parameters of type `Number`."
function parameters(m::T) where {T<:AbstractMaterial}
    Tuple([getfield(f, n) for n in fieldlabels(T) if fieldtype(T, n) isa Number])
end

label(m::AbstractMaterial) = m.label


#===========================#
# AbstractElement Materials #
#===========================#

include("./SVK.jl")

end # module