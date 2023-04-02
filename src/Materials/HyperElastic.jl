using .Materials: AbstractHyperElasticMaterial
using ..Utils: label, eye, _vogit
using Tensors: SymmetricTensor, hessian

import .Materials: density, parameters, cosserat, strain_energy

export HyperElastic

""" HyperElastic material struct.
### Fields:
- `params`         -- strain energy material parameters stored in a `Vector`.
- `Ψ(𝔼,params...)` -- strain energy function given a `Vector` of parameters `params` and 
                    the Green-Lagrange strain tensor `𝔼`.
- `ρ`              -- density (`nothing` for static cases).
- `label`          -- material label.

[See this ref.](https://en.wikipedia.org/wiki/Hyperelastic_material)
"""
struct HyperElastic{T<:Real,F<:Function,R<:Union{T,Nothing}} <: AbstractHyperElasticMaterial
    params::Vector{T}
    Ψ::F
    ρ::R
    label::Symbol
    function HyperElastic(params::Vector{T}, Ψ::F, ρ::R, label::L=:no_labelled_mat) where
    {T<:Real,F<:Function,R<:Union{Nothing,Real},L<:Union{Symbol,String}}
        return new{T,F,R}(params, Ψ, ρ, Symbol(label))
    end
end

"Constructor for an `HyperElastic` material with no density parameter `ρ`."
function HyperElastic(params::Vector{<:Real}, Ψ::Function, label::L=:no_labelled_mat) where {L<:Union{Symbol,String}}
    return HyperElastic(params, Ψ, nothing, label)
end


"Returns the strain energy function `Ψ` for a `HyperElastic` material `m`."
strain_energy(m::HyperElastic) = m.Ψ

"Returns the strain energy parameters `params` for a `HyperElastic` material `m`."
parameters(m::HyperElastic) = m.params

"Returns the Cosserat or Second-Piola Kirchoff stress tensor `𝕊` 
considering a `SVK` material `m` and the Lagrangian Green 
strain tensor `𝔼`.Also this function provides `∂𝕊∂𝔼` for the iterative method."
function cosserat(m::HyperElastic, 𝔼::AbstractMatrix)

    𝔼 = SymmetricTensor{2,3}(𝔼)

    # Closure strain energy function
    Ψ = E -> strain_energy(m)(E, parameters(m)...)

    ∂²Ψ∂E², 𝕊 = hessian(Ψ, 𝔼, :all)

    # Fill ∂S∂𝔼 with Belischko nomenclature
    ∂𝕊∂𝔼 = zeros(6, 6)
    indexes = [(1, 1), (2, 2), (3, 3), (2, 3), (1, 3), (1, 2)]

    row = 1
    for index in indexes
        i, j = index
        ∂𝕊∂𝔼[row, :] .= _vogit(∂²Ψ∂E²[:, :, i, j])
        row += 1
    end

    return 𝕊, ∂𝕊∂𝔼

end