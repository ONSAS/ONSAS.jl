using LinearAlgebra: tr
using SparseArrays: SparseMatrixCSC

using .Materials: AbstractMaterial
using ..Utils: label, eye

import .Materials: parameters, cosserat

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
function SVK(; λ::Real, G::Real, ρ::R=nothing, label::L=:no_labelled_mat) where
{R<:Union{Nothing,Real},L<:Union{Symbol,String}}

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


"Returns the Cosserat or Second-Piola Kirchoff tensor (𝕊) for a `Tetrahedron` element `t`
considering a `SVK` material `m` and the strain tensor `𝔼`."
function cosserat(m::SVK, 𝔼::AbstractMatrix, compute∂𝕊∂𝔼::Bool=true)

    λ, G = lame_parameters(m)
    𝕊 = λ * tr(𝔼) * eye(3) + 2 * G * 𝔼

    if compute∂𝕊∂𝔼
        ∂𝕊∂𝔼 = SparseMatrixCSC(zeros(6, 6))
        ∂𝕊∂𝔼[1:3, 1:3] = λ * ones(3, 3) + 2 * G * eye(3)
        ∂𝕊∂𝔼[4:6, 4:6] = G * eye(3)
        return 𝕊, ∂𝕊∂𝔼
    else
        return 𝕊
    end

end


