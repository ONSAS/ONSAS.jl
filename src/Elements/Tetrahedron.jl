using StaticArrays: SVector
using LinearAlgebra: det, diagm

using ..Materials: AbstractMaterial, SVK, cosserat
using ..Elements: AbstractElement, AbstractNode
using ..CrossSections: AbstractCrossSection, area
using ..Utils: eye, _vogit

import ..Elements: create_entity, internal_forces, local_dof_symbol, strain, stress

export Tetrahedron, volume

"""
A `Tetrahedron` represents a 3D volume element with four nodes.
### Fields:
- `nodes`    -- stores the tetrahedron nodes.
- `label`    -- stores the tetrahedron label.
"""
struct Tetrahedron{dim,T<:Real,N<:AbstractNode{dim,T}} <: AbstractElement{dim,T}
  nodes::SVector{4,N}
  label::Symbol
  function Tetrahedron(nodes::SVector{4,N}, label=:no_labelled_element) where
  {dim,T<:Real,N<:AbstractNode{dim,T}}
    @assert dim == 3 "Nodes of a tetrahedron element must be 3D."
    new{dim,T,N}(nodes, Symbol(label))
  end
end

Tetrahedron(n₁::N, n₂::N, n₃::N, n₄::N, label=:no_labelled_element) where {N<:AbstractNode} =
  Tetrahedron(SVector(n₁, n₂, n₃, n₄), Symbol(label))

"Constructor for a `Tetrahedron` element without nodes and a `label`. This function is used to create meshes via GMSH."
Tetrahedron(label::L=:no_labelled_face) where {L<:Union{String,Symbol}} = 
Tetrahedron(SVector(Node(0,0,0),Node(0,0,0),Node(0,0,0),Node(0,0,0)), Symbol(label))

"Returns a `Tetrahedron` given an empty `Tetrahedron` `t` and a `Vector` of `Node`s `vn`."
create_entity(t::Tetrahedron, vn::AbstractVector{<:AbstractNode}) = Tetrahedron(vn[1], vn[2], vn[3], vn[4], label(t))


"Returns the `Tetrahedron` `t` volume in the reference configuration."
function volume(t::Tetrahedron)
  d = _shape_functions_derivatives(t)
  coords = _coordinates_matrix(t)
  J = _jacobian_mat(coords, d)
  vol = _volume(J)
  return vol
end


"Returns the local dof symbol of a `Tetrahedron` element."
local_dof_symbol(::Tetrahedron) = [:u]

"Returns the shape functions derivatives of a `Tetrahedron` element."
function _shape_functions_derivatives(::Tetrahedron, order =1)
  d = zeros(3, 4)
  if order == 1
    d[1, 1] = 1
    d[1:3, 2] = [-1, -1, -1]
    d[3, 3] = 1
    d[2, 4] = 1
  end
  return d
end

"Returns the reshaped coordinates `elem_coords` of the tetrahedron element into a 4x3 matrix."
_coordinates_matrix(t::Tetrahedron) = reduce(hcat,coordinates(t))

"Computes Jacobian matrix"
_jacobian_mat(tetrahedron_coords_matrix::AbstractMatrix, derivatives::AbstractMatrix) = tetrahedron_coords_matrix * derivatives'

"Computes volume element of a tetrahedron given J = det(𝔽)."
function _volume(jacobian_mat::AbstractMatrix)
  volume = det(jacobian_mat) / 6.0
  @assert volume > 0 throw(ArgumentError("Element with negative volume, check connectivity."))
  return volume
end


function _B_mat(deriv::AbstractMatrix , 𝔽::AbstractMatrix)

  B = zeros(6, 12) 

  B[1:3, :] = [diagm(deriv[:,1])*𝔽' diagm(deriv[:,2])*𝔽' diagm(deriv[:,3])*𝔽' diagm(deriv[:,4])*𝔽']

  for k in 1:4
      B[4:6 , (k-1)*3 .+ (1:3)] = [ deriv[2,k] * 𝔽[:,3]' + deriv[3,k] * 𝔽[:,2]'
                                    deriv[1,k] * 𝔽[:,3]' + deriv[3,k] * 𝔽[:,1]'
                                    deriv[1,k] * 𝔽[:,2]' + deriv[2,k] * 𝔽[:,1]' ] 
  end
  return B
end

"Returns the internal force of a `Tetrahedron` element `t` doted with an `AbstractMaterial` `m` +
and a an element displacement vector `u_e`."
function internal_forces(f)

  d = _shape_functions_derivatives(t)

  coords = _coordinates_matrix(t)

  disps = reshape(u_e, 3, 4)

  J = _jacobian_mat(coords, d)

  vol = _volume(J)

  funder = inv(J)' * d
  
  # ∂u∂X in global coordinats 
  ℍ = disps * funder'

  # Deformation gradient 
  𝔽 = ℍ + eye(3)

  # Cauchy strain tensor
  ℂ = 𝔽' * 𝔽 

  # Green-Lagrange strain  
  𝔼 = 0.5 * (ℍ + ℍ' + ℍ' * ℍ)

  𝕊, ∂𝕊∂𝔼 = cosserat(m, 𝔼)

  B = _B_mat(funder, 𝔽)

  𝕊_vogit = _vogit(𝕊)

  fᵢₙₜ_e = B' * 𝕊_vogit * vol
  
  # Material stiffness
  Kₘ = B' * ∂𝕊∂𝔼 * B* vol

  # Geometric stiffness
  aux = funder' * 𝕊 * funder  * vol 

  Kᵧ = zeros(12,12) 

  for i in 1:4
    for j in 1:4
      Kᵧ[(i-1)*3+1 , (j-1) * 3 + 1] = aux[i,j]
      Kᵧ[(i-1)*3+2 , (j-1) * 3 + 2] = aux[i,j]
      Kᵧ[(i-1)*3+3 , (j-1) * 3 + 3] = aux[i,j]
    end
  end

  # Stifness matrix
  Kᵢₙₜ_e = Kₘ + Kᵧ

  # Piola stress
  ℙ = 𝔽 * 𝕊

  # Cuachy stress
  # σ_e = ℙ
  
  return fᵢₙₜ_e, Kᵢₙₜ_e, ℙ, ℂ

end




