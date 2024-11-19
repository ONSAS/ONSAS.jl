"Module defining a generic hyper elastic material model form the strain energy function."
module HyperElasticMaterial

using Tensors, Reexport
using ..HyperElasticMaterials
using ..Utils

@reexport import ..Materials: parameters
@reexport import ..HyperElasticMaterials: cosserat_stress!, strain_energy

export HyperElastic

"""
Material with hyperelastic properties.

For context see the wikipedia article on [Hyperelastic_material](https://en.wikipedia.org/wiki/Hyperelastic_material).
"""
struct HyperElastic{T <: Real, F <: Function} <: AbstractHyperElasticMaterial
    "Strain energy material parameters."
    params::Vector{T}
    "Strain energy function given `params` and the Green-Lagrange strain tensor `𝔼`."
    Ψ::F
    "Density (`nothing` for static cases)."
    ρ::Density
    "Material label."
    label::Label
    function HyperElastic(params::Vector{T}, Ψ::F, ρ::Density,
            label::Label = NO_LABEL) where {T <: Real, F <: Function}
        new{T, F}(params, Ψ, ρ, Symbol(label))
    end
end

"Constructor for `HyperElastic` material with no density."
function HyperElastic(
        params::Vector{T}, Ψ::F, label::Label = NO_LABEL) where {T <: Real, F <: Function}
    HyperElastic(params, Ψ, nothing, label)
end

"Return the strain energy function `Ψ` for a `HyperElastic` material `m`."
strain_energy(m::HyperElastic) = m.Ψ

"Return the strain energy parameters `params` for a `HyperElastic` material `m`."
parameters(m::HyperElastic) = m.params

"Return the Cosserat or Second-Piola Kirchoff stress tensor `𝕊`
considering a `SVK` material `m` and the Lagrangian Green
strain tensor `𝔼`.Also this function provides `∂𝕊∂𝔼` for the iterative method."
function cosserat_stress!(S::AbstractMatrix{<:Real}, ∂S∂E::AbstractMatrix{<:Real},
        m::HyperElastic, E::AbstractMatrix)
    # Transform 𝔼 to a Tenor
    𝔼 = SymmetricTensor{2, 3}(E)
    # Closure strain energy function
    Ψ = E -> strain_energy(m)(E, parameters(m)...)

    ∂²Ψ∂E², 𝕊 = hessian(Ψ, 𝔼, :all)
    # Fill symmetric matrix
    fill_symmetric_matrix!(S, 𝕊[1, 1], 𝕊[2, 2], 𝕊[3, 3], 𝕊[2, 3], 𝕊[1, 3], 𝕊[1, 2])

    # Fill ∂S∂𝔼 with Belischko nomenclature
    row = 1
    for (i, j) in INDEXES_TO_VOIGT
        ∂S∂E[row, :] .= voigt(∂²Ψ∂E²[:, :, i, j])
        row += 1
    end
    S, ∂S∂E
end

end
