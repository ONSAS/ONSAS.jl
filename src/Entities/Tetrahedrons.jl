"Module defining tetrahedron elements."
module Tetrahedrons

using StaticArrays, LinearAlgebra, LazySets, Reexport

using ..Utils
using ..Nodes
using ..Entities
using ..Materials
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

"Contains the cache to compute the element internal forces and stiffness matrix."
struct TetrahedronCache{T,ST<:Symmetric{T}} <: AbstractElementCache
    "Internal forces."
    fint::Vector{T}
    "Stiffness matrix."
    Ks::ST
    "Piola stress."
    σ::ST
    "Cauchy-Green strain."
    ε::ST
    "Deformation gradient."
    F::Matrix{T}
    H::Matrix{T}
    "Reference coordinates."
    X::Matrix{T}
    "Jacobian matrix."
    J::Matrix{T}
    "Shape functions derivatives."
    funder::Matrix{T}
    "B matrix."
    B::Matrix{T}
    "Auxiliary matrix for computing Geometric Stiffness (only AbstractHyperElasticMaterial)."
    aux_geometric_Ks::Matrix{T}
    "Lagrange Green Strain"
    E::ST
    function TetrahedronCache()
        fint = zeros(12)
        Ks = Symmetric(zeros(12, 12))
        σ = Symmetric(zeros(3, 3))
        ε = Symmetric(zeros(3, 3))
        F = zeros(3, 3)
        H = zeros(3, 3)
        X = zeros(3, 4)
        J = zeros(3, 3)
        funder = zeros(3, 4)
        B = zeros(6, 12)
        aux_geometric_Ks = zeros(4, 4)
        E = Symmetric(zeros(3, 3))
        new{Float64,Symmetric{Float64}}(fint, Ks, σ, ε, F, H, X, J,
                                        funder, B, aux_geometric_Ks, E)
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
function _jacobian_mat(tetrahedron_coords_matrix::AbstractMatrix, derivatives::Matrix)
    tetrahedron_coords_matrix * derivatives'
end

"Computes volume element of a tetrahedron given J = det(F)."
function _volume(jacobian_mat::Matrix)
    volume = det(jacobian_mat) / 6.0
    @assert volume > 0 throw(ArgumentError("Element with negative volume, check connectivity."))
    volume
end

function _B_mat!(B::Matrix, deriv::Matrix, F::Matrix)
    B[1:3, :] = [diagm(deriv[:, 1]) * F' diagm(deriv[:, 2]) *
                                         F' diagm(deriv[:, 3]) * F' diagm(deriv[:, 4]) * F']

    for k in 1:4
        B[4:6, (k - 1) * 3 .+ (1:3)] = [deriv[2, k] * F[:, 3]' + deriv[3, k] * F[:, 2]'
                                        deriv[1, k] * F[:, 3]' + deriv[3, k] * F[:, 1]'
                                        deriv[1, k] * F[:, 2]' + deriv[2, k] * F[:, 1]']
    end
    B
end

"Return the internal forces of a `Tetrahedron` element `t` doted with and `AbstractMaterial` `m` and
an element displacement vector `u_e`. "
function internal_forces(m::AbstractMaterial, t::Tetrahedron, u_e::AbstractVector)
    internal_forces(m, t, u_e, TetrahedronCache())
end

"Return the geometric stiffness."
function geometric_stiffness!(Ks::Symmetric, aux_geometric_Ks::Matrix,
                              𝕊::AbstractMatrix, funder::Matrix{<:Real}, vol::Real)
    aux_geometric_Ks .= funder' * 𝕊 * funder * vol
    for i in 1:4
        for j in i:4
            # Diagonal elements can be set directly
            if i == j
                Ks[(i - 1) * 3 + 1, (j - 1) * 3 + 1] = aux_geometric_Ks[i, j]
                Ks[(i - 1) * 3 + 2, (j - 1) * 3 + 2] = aux_geometric_Ks[i, j]
                Ks[(i - 1) * 3 + 3, (j - 1) * 3 + 3] = aux_geometric_Ks[i, j]
                # Off-diagonal elements must ensure symmetry
            else
                Ks.data[(i - 1) * 3 + 1, (j - 1) * 3 + 1] = aux_geometric_Ks[i, j]
                Ks.data[(j - 1) * 3 + 1, (i - 1) * 3 + 1] = aux_geometric_Ks[i, j]
                Ks.data[(i - 1) * 3 + 2, (j - 1) * 3 + 2] = aux_geometric_Ks[i, j]
                Ks.data[(j - 1) * 3 + 2, (i - 1) * 3 + 2] = aux_geometric_Ks[i, j]
                Ks.data[(i - 1) * 3 + 3, (j - 1) * 3 + 3] = aux_geometric_Ks[i, j]
                Ks.data[(j - 1) * 3 + 3, (i - 1) * 3 + 3] = aux_geometric_Ks[i, j]
            end
        end
    end
end

"Return the internal force of a `Tetrahedron` element `t` doted with an `AbstractHyperElasticMaterial` `m` +
and a an element displacement vector `u_e`. This function modifies the cache to avoid memory allocations."
function internal_forces(m::AbstractHyperElasticMaterial, t::Tetrahedron, u_e::AbstractVector,
                         cache::TetrahedronCache)
    (; fint, Ks, σ, ε, F, H, X, J, funder, B, aux_geometric_Ks, E) = cache

    # Kinematics
    U = reshape(u_e, 3, 4)
    ∂X∂ζ = _shape_functions_derivatives(t)
    X .= _coordinates_matrix(t)
    J .= _jacobian_mat(X, ∂X∂ζ)
    vol = _volume(J)
    funder .= inv(J)' * ∂X∂ζ
    H .= U * funder'
    F .= H + eye(3)

    E .= Symmetric(0.5 * (H + H' + H' * H))
    B = _B_mat!(B, funder, F)

    # Stresses
    𝕊, ∂𝕊∂𝔼 = cosserat_stress(m, E)
    𝕊_voigt = voigt(𝕊)
    fint .= B' * 𝕊_voigt * vol

    # Material stiffness
    Km = Symmetric(B' * ∂𝕊∂𝔼 * B * vol)

    # Geometric stiffness
    Ks .= 0.0
    geometric_stiffness!(Ks, aux_geometric_Ks, 𝕊, funder, vol)
    Ks .= Km + Ks

    # Piola stress
    σ .= Symmetric(F * 𝕊)

    # Right hand Cauchy strain tensor
    ε .= Symmetric(F' * F)

    fint, Ks, σ, ε
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
function internal_forces(m::IsotropicLinearElastic, t::Tetrahedron, u_e::AbstractVector,
                         cache::TetrahedronCache)
    (; fint, Ks, σ, ε, F, H, X, J, funder, B) = cache

    # Kinematics
    ∂X∂ζ = _shape_functions_derivatives(t)
    X .= _coordinates_matrix(t)
    J .= _jacobian_mat(X, ∂X∂ζ)
    vol = _volume(J)
    funder .= inv(J)' * ∂X∂ζ
    U = reshape(u_e, 3, 4)
    H .= U * funder'
    ε .= Symmetric(0.5 * (H + H'))
    F .= eye(3)
    _B_mat!(B, funder, F)

    # Stresses
    σaux, ∂𝕊∂𝔼 = cauchy_stress(m, ε)
    σ .= Symmetric(σaux)

    # Stiffness matrix
    Ks .= Symmetric(B' * ∂𝕊∂𝔼 * B * vol)

    fint .= Ks * u_e

    fint, Ks, σ, ε
end

"Shape function derivatives."
const ∂X∂ζ_1 = [1.0  -1.0  0.0  0.0
                0.0  -1.0  0.0  1.0
                0.0  -1.0  1.0  0.0]

"Return the shape functions derivatives of a `Tetrahedron` element."
function _shape_functions_derivatives(::Tetrahedron, order::Int=1)
    ∂X∂ζ = if order == 1
        ∂X∂ζ_1
    end
end

"Indices for computing the minors of the interpolation matrix, implemented as a hash table."
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
