"Module defining tetrahedron elements."
module Tetrahedrons

using StaticArrays, LinearAlgebra, LazySets, Reexport

using ..Utils
using ..Nodes
using ..Entities
using ..IsotropicLinearElasticMaterial
using ..HyperElasticMaterials

@reexport import ..Entities: create_entity, internal_forces, local_dof_symbol, strain, stress,
                             weights, volume, elements_cache

export Tetrahedron, reference_coordinates

"""
A `Tetrahedron` represents a 3D volume element with four nodes.

See [[Belytschko]](@ref) and [[Gurtin]](@ref) for more details.
"""
struct Tetrahedron{dim,T<:Real,N<:AbstractNode{dim,T},VN<:AbstractVector{N}} <:
       AbstractElement{dim,T}
    "Tetrahedron nodes."
    nodes::VN
    "Tetrahedron label."
    label::Label
    function Tetrahedron(nodes::VN,
                         label::Label=NO_LABEL) where
             {dim,T<:Real,N<:AbstractNode{dim,T},VN<:AbstractVector{N}}
        @assert dim == 3 "Nodes of a tetrahedron element must be 3D."
        new{dim,T,N,VN}(nodes, Symbol(label))
    end
end
function Tetrahedron(n₁::N, n₂::N, n₃::N, n₄::N, label::Label=NO_LABEL) where {N<:AbstractNode}
    Tetrahedron(SVector(n₁, n₂, n₃, n₄), label)
end
function Tetrahedron(nodes::AbstractVector{N}, label::Label=NO_LABEL) where {N<:AbstractNode}
    Tetrahedron(SVector(nodes...), label)
end
"Constructor for a `Tetrahedron` element without nodes and a `label`. This function is used to create meshes via GMSH."
function Tetrahedron(label::Label=NO_LABEL)
    Tetrahedron(SVector(Node(0, 0, 0), Node(0, 0, 0), Node(0, 0, 0), Node(0, 0, 0)), label)
end

"Return a `Tetrahedron` given an empty `Tetrahedron` `t` and a `Vector` of `Node`s `vn`."
function create_entity(t::Tetrahedron, vn::AbstractVector{<:AbstractNode})
    Tetrahedron(vn, label(t))
end

struct TetrahedronCache{T} <: AbstractElementCache
    fint::Vector{T}
    # TODO Use Symmetric K.
    Ks::Matrix{T}
    σ::Matrix{T}
    ε::Matrix{T}
    F::Matrix{T}
    H::Matrix{T}
    function TetrahedronCache()
        fint = zeros(12)
        Ks = zeros(12, 12)
        σ = zeros(3, 3)
        ε = zeros(3, 3)
        F = zeros(3, 3)
        H = zeros(3, 3)
        new{Float64}(fint, Ks, σ, ε, F, H)
    end
end

"Return an empty Tetrahedron cache."
elements_cache(::Type{Tetrahedron}) = TetrahedronCache()

"Return the `Tetrahedron` `t` volume in the reference configuration."
function volume(t::Tetrahedron)
    ∂X∂ζ = _shape_functions_derivatives(t)
    coords = _coordinates_matrix(t)
    J = _jacobian_mat(coords, ∂X∂ζ)
    vol = _volume(J)
end

"Return the local dof symbol of a `Tetrahedron` element."
local_dof_symbol(::Tetrahedron) = [:u]

"Return the reshaped coordinates `elem_coords` of the tetrahedron element into a 4x3 matrix."
_coordinates_matrix(t::Tetrahedron) = reduce(hcat, coordinates(t))

"Computes Jacobian matrix"
function _jacobian_mat(tetrahedron_coords_matrix::AbstractMatrix, derivatives::AbstractMatrix)
    tetrahedron_coords_matrix * derivatives'
end

"Computes volume element of a tetrahedron given J = det(𝔽)."
function _volume(jacobian_mat::AbstractMatrix)
    volume = det(jacobian_mat) / 6.0
    @assert volume > 0 throw(ArgumentError("Element with negative volume, check connectivity."))
    volume
end

function _B_mat(deriv::AbstractMatrix, 𝔽::AbstractMatrix)
    B = zeros(6, 12)

    B[1:3, :] = [diagm(deriv[:, 1]) * 𝔽' diagm(deriv[:, 2]) * 𝔽' diagm(deriv[:, 3]) * 𝔽' diagm(deriv[:,
                                                                                                     4]) *
                                                                                         𝔽']

    for k in 1:4
        B[4:6, (k - 1) * 3 .+ (1:3)] = [deriv[2, k] * 𝔽[:, 3]' + deriv[3, k] * 𝔽[:, 2]'
                                        deriv[1, k] * 𝔽[:, 3]' + deriv[3, k] * 𝔽[:, 1]'
                                        deriv[1, k] * 𝔽[:, 2]' + deriv[2, k] * 𝔽[:, 1]']
    end
    B
end

function internal_forces(m::AbstractHyperElasticMaterial, t::Tetrahedron, u_e::AbstractVector)
    internal_forces(m, t, u_e, TetrahedronCache())
end

"Return the internal force of a `Tetrahedron` element `t` doted with an `AbstractHyperElasticMaterial` `m` +
and a an element displacement vector `u_e`."
function internal_forces(m::AbstractHyperElasticMaterial, t::Tetrahedron, u_e::AbstractVector,
                         cache::TetrahedronCache)
    (; fint, Ks, σ, ε, F, H) = cache

    ∂X∂ζ = _shape_functions_derivatives(t)
    X = _coordinates_matrix(t)
    U = reshape(u_e, 3, 4)
    J = _jacobian_mat(X, ∂X∂ζ)
    vol = _volume(J)

    # The deformation gradient F can be obtained by integrating f under over time ∂F/∂t.
    funder = inv(J)' * ∂X∂ζ

    # ∇u in global coordinats
    H .= U * funder'

    # Deformation gradient
    F .= H + eye(3)

    # Green-Lagrange strain
    𝔼 = Symmetric(0.5 * (H + H' + H' * H))
    𝕊, ∂𝕊∂𝔼 = cosserat_stress(m, 𝔼)
    B = _B_mat(funder, F)
    𝕊_voigt = voigt(𝕊)
    fint .= B' * 𝕊_voigt * vol

    # Material stiffness
    Kₘ = Symmetric(B' * ∂𝕊∂𝔼 * B * vol)
    Ks .= 0.0

    # Geometric stiffness
    aux = funder' * 𝕊 * funder * vol
    for i in 1:4
        for j in 1:4
            Ks[(i - 1) * 3 + 1, (j - 1) * 3 + 1] = aux[i, j]
            Ks[(i - 1) * 3 + 2, (j - 1) * 3 + 2] = aux[i, j]
            Ks[(i - 1) * 3 + 3, (j - 1) * 3 + 3] = aux[i, j]
        end
    end

    # Stifness matrix
    Ks .= Kₘ + Ks

    # Compute stress and strian just for post-process
    # Piola stress
    σ .= Symmetric(F * 𝕊)

    # Cauchy strain tensor
    ε .= Symmetric(F' * F)

    fint, Ks, σ, ε
end

# TODO Implement method.
function internal_forces(m::IsotropicLinearElastic, t::Tetrahedron, u_e::AbstractVector,
                         ::TetrahedronCache)
    internal_forces(m, t, u_e)
end

"
Return the internal force of a `Tetrahedron` element `t` doted with an `LinearIsotropicMaterial` `m`.
## Arguments
- `material`: `IsotropicLinearElastic` type, the linear elastic material of the tetrahedron element.
- `element`: `Tetrahedron` type, the tetrahedron element for which internal forces are to be computed.
- `displacements`: `AbstractVector` type, the nodal displacements of the element.

## Return
A 4-tuple containing:
- `forces`: `Symmetric` type, the internal forces of the tetrahedron element.
- `stiffness`: `Symmetric` type, the stiffness matrix of the tetrahedron element.
- `stress`: `Symmetric` type, the Cauchy stress tensor of the tetrahedron element.
- `strain`: `Symmetric` type, the strain tensor of the tetrahedron element.
"
function internal_forces(m::IsotropicLinearElastic, t::Tetrahedron, u_e::AbstractVector)
    ∂X∂ζ = _shape_functions_derivatives(t)

    X = _coordinates_matrix(t)

    J = _jacobian_mat(X, ∂X∂ζ)

    vol = _volume(J)

    funder = inv(J)' * ∂X∂ζ

    # ∇u = H in global coordinats
    U = reshape(u_e, 3, 4)
    H = U * funder'

    ϵ = Symmetric(0.5 * (H + H'))
    𝔽 = eye(3)

    B = _B_mat(funder, 𝔽)

    σ, ∂σ∂ϵ = cauchy_stress(m, ϵ)

    Kᵢₙₜ_e = Symmetric(B' * ∂σ∂ϵ * B * vol)

    fᵢₙₜ_e = Kᵢₙₜ_e * u_e

    fᵢₙₜ_e, Kᵢₙₜ_e, σ, ϵ
end

"Return the shape functions derivatives of a `Tetrahedron` element."
function _shape_functions_derivatives(::Tetrahedron, order=1)
    d = zeros(3, 4)
    if order == 1
        d[1, 1] = 1
        d[1:3, 2] = [-1, -1, -1]
        d[3, 3] = 1
        d[2, 4] = 1
    end
    d
end

"""
Indices for computing the minors of the interpolation matrix, implemented as a hash table.
"""
const MINOR_INDICES = [([2, 3, 4], [2, 3, 4])    ([2, 3, 4], [1, 3, 4])    ([2, 3, 4], [1, 2, 4])    ([2, 3, 4], [1, 2, 3])
                       ([1, 3, 4], [2, 3, 4])    ([1, 3, 4], [1, 3, 4])    ([1, 3, 4], [1, 2, 4])    ([1, 3, 4], [1, 2, 3])
                       ([1, 2, 4], [2, 3, 4])    ([1, 2, 4], [1, 3, 4])    ([1, 2, 4], [1, 2, 4])    ([1, 2, 4], [1, 2, 3])
                       ([1, 2, 3], [2, 3, 4])    ([1, 2, 3], [1, 3, 4])    ([1, 2, 3], [1, 2, 4])    ([1, 2, 3], [1, 2, 3])]

"Return the interpolation matrix `𝑀` for a `Tetrahedron` element `t`."
function interpolation_matrix(t::Tetrahedron{3,T}) where {T<:Real}
    # Node coordinates matrix 𝐴.
    𝐴 = MMatrix{4,4,T}(undef)

    @inbounds for (node_index, node) in enumerate(nodes(t))
        𝐴[node_index, 1] = one(T)
        𝐴[node_index, 2:4] = coordinates(node)
    end

    # 𝑀 matrix.
    𝑀 = MMatrix{4,4,T}(undef)
    V = det(𝐴)

    # Compute minors.
    @inbounds for J in 1:4, I in 1:4
        ridx, cidx = MINOR_INDICES[I, J]
        Â = det(𝐴[ridx, cidx])
        𝑀[I, J] = Â / V * (-1)^(I + J)
    end
    𝑀
end

"Return the interpolation weights of a point `p` in a `Tetrahedron` element `t`."
function weights(t::Tetrahedron{3}, p::Point{3})
    interpolation_matrix(t) * Point(1, p...)
end

function Base.convert(::Type{LazySets.Tetrahedron}, t::Tetrahedron)
    LazySets.Tetrahedron(nodes(t))
end

"Checks if a point `p` is inside a `Tetrahedron` element `t`."
Base.:∈(p::Point, t::Tetrahedron) = p ∈ convert(LazySets.Tetrahedron, t)

end
