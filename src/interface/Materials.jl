#######################
# Materials interface #
#######################

"""
Module defining the materials implemented.
Each material consists of a data type with one or more parameters into its fields.
"""
module Materials

using Reexport: @reexport

@reexport import ..Utils: ScalarWrapper, label, set_label!

export AbstractMaterial, parameters
export SVK, lame_parameters


""" Abstract supertype for all material models.

An `AbstractMaterial` object facilitates the process of defining new material models. 
Different material models leads to different constitutive laws, internal forces and stiffness matrices.


**Common methods:**

* [`parameters`](@ref)
* [`label`](@ref)
* [`set_label!`](@ref)
"""
abstract type AbstractMaterial end

"Returns the parameters of type `Number`."
function parameters(m::T) where {T<:AbstractMaterial}
    Tuple([getfield(f, n) for n in fieldlabels(T) if fieldtype(T, n) isa Number])
end

label(m::AbstractMaterial) = m.label[]
set_label!(m::AbstractMaterial, label::Symbol) = m.label[] = label
set_label!(m::AbstractMaterial, label::String) = set_label!(m, Symbol(label))

const DEFAULT_LABEL = :label_no_assignned

#############################
# Materials implementations #
#############################

""" SVK material struct.
### Fields:
- `E` -- Elasticity modulus.
- `ν` -- Poisson's ratio.
- `ρ` -- Density (`nothing` for static cases).
- `label` -- Label of the material

"""
struct SVK <: AbstractMaterial
    E::Number
    ν::Number
    ρ::Union{<:Number,Nothing}
    label::ScalarWrapper{Symbol}
    function SVK(E, ν, ρ=nothing, label=DEFAULT_LABEL)

        return new(E, ν, ρ, ScalarWrapper(Symbol(label)))
    end
end

"Constructor with lamé parameters λ and G"
function SVK(; λ, G, ρ, label)

    # Compute E and ν given Lamé parameters λ and μ (μ = G)
    E = G * (3λ + 2G) / (λ + G)
    ν = λ / (2(λ + G))

    return SVK(E, ν, ρ, Symbol(label))
end


"Returns SVK parameters tuple."
parameters(m::SVK) = (m.E, m.ν)

" Convert svk material parameters to lamé nomenclature"
function lame_parameters(svk::SVK)

    E = svk.E
    ν = svk.ν

    # Compute Lamé parameters λ and G
    G = E / (2(1 + ν))
    λ = E * ν / ((1 + ν) * (1 - 2 * ν))

    return λ, G
end

end # module