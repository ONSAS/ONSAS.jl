using .Materials: AbstractMaterial
using ..Utils: label

import .Materials: parameters

export SVK, lame_parameters

""" SVK material struct.
### Fields:
- `E` -- Elasticity modulus.
- `ν` -- Poisson's ratio.
- `ρ` -- Density (`nothing` for static cases).
- `label` -- Label of the material
"""
struct SVK{T<:Real} <: AbstractMaterial
    E::T
    ν::T
    ρ::Union{T,Nothing}
    label::Symbol
    function SVK(E::T, ν::T, ρ::R, label::L=:no_labelled_mat) where
    {T<:Real,R<:Union{Nothing,Real},L<:Union{Symbol,String}}
        return new{T}(E, ν, ρ, Symbol(label))
    end
end

"Material `SVK` constructor with no density parameter `ρ`"
function SVK(E::Real, ν::Real, label::L=:no_labelled_mat) where {L<:Union{Symbol,String}}
    return SVK(E, ν, nothing, label)
end

"Material `SVK` constructor lamé parameters `λ` and `G`"
function SVK(; λ::Real, G::Real, ρ::Real, label)

    # Compute E and ν given Lamé parameters λ and μ (μ = G)
    E = G * (3λ + 2G) / (λ + G)
    ν = λ / (2(λ + G))

    return SVK(E, ν, ρ, Symbol(label))
end

"Returns material `m` SVK parameters into a tuple."
parameters(m::SVK) = (m.E, m.ν)

"Returns lamé parameters `λ` and `G` from a `SVK`` material `m`."
function lame_parameters(svk::SVK)

    E = svk.E
    ν = svk.ν

    # Compute Lamé parameters λ and G
    G = E / (2(1 + ν))
    λ = E * ν / ((1 + ν) * (1 - 2 * ν))

    return λ, G
end
