# Implement a linear elastic material
using .Materials: AbstractMaterial

import .Materials: lame_parameters, shear_modulus, poisson_ratio, elasticity_modulus, bulk_modulus

export IsotropicLinearElastic

""" IsotropicLinearElastic material struct.
### Fields:
- `E`         -- Elasticity modulus.
- `ν`         -- Poisson's ratio.
- `ρ`         -- density (`nothing` for static cases).
- `label`     -- material label.

[See this ref.](https://en.wikipedia.org/wiki/Hyperelastic_material)
"""
struct IsotropicLinearElastic{ET<:Real,NT<:Real,RT<:Union{ET,Nothing}} <: AbstractMaterial
    E::ET
    ν::NT
    ρ::RT
    label::Symbol
    function IsotropicLinearElastic(E::ET, ν::NT, ρ::RT, label::L=:no_labelled_mat) where
    {ET<:Real,NT<:Real,RT<:Union{Nothing,Real},L<:Union{Symbol,String}}
        return new{ET,NT,RT}(E, ν, ρ, Symbol(label))
    end
end
"Material `IsotropicLinearElastic` constructor with no density parameter `ρ`."
IsotropicLinearElastic(E::Real, ν::Real, label::L=:no_labelled_mat) where {L<:Union{Symbol,String}} =
    IsotropicLinearElastic(E, ν, nothing, label)

"Material `IsotropicLinearElastic` constructor with Lamé parameters `λ`, `G` and density `ρ`. 
See [this ref](https://en.wikipedia.org/wiki/Lam%C3%A9_parameters)."
function IsotropicLinearElastic(; λ::Real, G::Real, ρ::R=nothing, label::L=:no_labelled_mat) where
{R<:Union{Nothing,Real},L<:Union{Symbol,String}}

    E = G * (3λ + 2G) / (λ + G)
    ν = λ / (2 * (λ + G))

    IsotropicLinearElastic(E, ν, ρ, Symbol(label))
end

"Returns the Elasticity modulus `E` from a `IsotropicLinearElastic` material `m`."
elasticity_modulus(m::IsotropicLinearElastic) = m.E

"Returns the Poisson's ratio `ν` from a `IsotropicLinearElastic` material `m`."
poisson_ratio(m::IsotropicLinearElastic) = m.ν

"Returns the shear modulus `G` from a `IsotropicLinearElastic` material `m`."
shear_modulus(m::IsotropicLinearElastic) = elasticity_modulus(m) / (2 * (1 + poisson_ratio(m)))

"Returns the bulk modulus `K` from a `IsotropicLinearElastic` material `m`."
bulk_modulus(m::IsotropicLinearElastic) = elasticity_modulus(m) / (3 * (1 - 2 * poisson_ratio(m)))
