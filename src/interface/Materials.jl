#######################
# Materials interface #
#######################

"""
Module defining the materials implemented.
Each material consists of a data type with one or more parameters into its fields.
"""
module Materials

export AbstractMaterial, model, label, parameters
export SVK, lame_parameters

""" Abstract supertype for a material.
The following methods are provided by the interface:
- `model(m)`      -- return a string with the material model (defaults to the material's type label)
- `parameters(m)` -- return a tuple with the material parameters
- `label(m)`      -- return material label
"""
abstract type AbstractMaterial end

"Returns the material model of a `T`` material type. "
model(::Type{T}) where {T<:AbstractMaterial} = string(T)

function parameters(m::T) where {T<:AbstractMaterial}
    Tuple([getfield(f, n) for n in fieldlabels(T) if fieldtype(T, n) isa Number])
end

"Returns the material with label `m`."
label(m::AbstractMaterial) = "no label is implemented, please overload this method."

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
    label::Symbol
    function SVK(E, ν, ρ=nothing, label=DEFAULT_LABEL)

        return new(E, ν, ρ, Symbol(label))
    end
end

"Constructor with lamé parameters λ and G"
function SVK(; λ, G, ρ, label)

    # Compute E and ν given Lamé parameters λ and μ (μ = G)
    E = G * (3λ + 2G) / (λ + G)
    ν = λ / (2(λ + G))

    return SVK(E, ν, ρ, Symbol(label))
end


"Returns SVK material model label."
model(::SVK) = "SVK"

"Returns SVK parameters tuple."
parameters(m::SVK) = (m.E, m.ν)

"Returns SVK label"
label(m::SVK) = string(m.label)

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