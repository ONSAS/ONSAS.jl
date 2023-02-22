using ..Materials: AbstractMaterial, SVK
using ..Utils: ScalarWrapper, eye
using StaticArrays: @SMatrix, @SVector
export Tetrahedron

"""
A `Tetrahedron` represents a 3D volume element.
### Fields:
- `nodes`          -- stores truss nodes.
- `material`       -- stores truss material.
- `label`          -- stores the truss label.
"""
struct Tetrahedron{dim,M} <: AbstractElement{dim,M}
  nodes::Vector{<:AbstractNode{dim}}
  material::M
  label::ScalarWrapper{Symbol}
end

#TODO: Add tetrahedron element order as a parametric type
"Derivatives of the linear shape functions (Order 1)"
function _shape_functions_derivatives()
  d = zeros(3, 4)
  d[1, 1] = 1
  d[1:3, 2] = [-1, -1, -1]
  d[3, 3] = 1
  d[2, 4] = 1
  return d
end

"Returns the reshaped coordinates of the tetrahedron element"
_tetra_coords_mat(elem_coords) = reshape(transpose(elem_coords), 3, 4)

"Computes Jacobian matrix"
_jacobian_mat(tetrahedron_coords_matrix::AbstractMatrix, derivatives::AbstractMatrix) = tetrahedron_coords_matrix * derivatives'

"Computes volume o"
function _volume(jacobian_mat::AbstractMatrix)
  volume = det(jacobian_mat) / 6.0
  volume < 0 && throw(ArgumentError("Element with negative volume, check connectivity."))
  return volume
end


"Computes tetrahedron stiffness matrix."
function _computes_cosserat_and_derivatives(t::Tetrahedron{3}, u_e::AbstractVector)

  d = _shape_functions_derivatives()

  tetra_coords = _tetra_coords_mat(coordinates(t))

  jacobian_mat = _jacobian_mat(tetra_coords, d)

  vol = _volume(jacobian_mat)

  funder = inv(jacobian_mat)' * d

  u_tetra_mat = reshape(u_e, 3, 4)

  H = u_tetra_mat * funder

  𝔽 = H + eye(e)

  𝔼 = 0.5 * (H + H' + H' * H)

  𝕊, ∂𝕊∂𝔼 = _cosserat_tensor(material(m), 𝔼)

  return 𝕊, ∂𝕊∂𝔼
end

"Computes tetrahedron internal force"
function internal_force(t::Tetrahedron{3}, u_e::AbstractVector)

  fₑ = @SVector zeros(12)

  d = _shape_functions_derivatives()

  tetra_coords = _tetra_coords_mat(coordinates(t))

  jacobian_mat = _jacobian_mat(tetra_coords, d)

  vol = _volume(jacobian_mat)

  funder = inv(jacobian_mat)' * d

  u_tetra_mat = reshape(u_e, 3, 4)

  H = u_tetra_mat * funder

  𝔽 = H + eye(e)

  𝔼 = 0.5 * (H + H' + H' * H)

  𝕊, ∂𝕊∂𝔼 = _cosserat_tensor(material(m), 𝔼)

  B = _compute_B_mat(funder, 𝔽)

  𝕊_vogit = matrix2vogit(𝕊)

  # Internal force
  fₑ = SVector(B' * 𝕊_vogit * vol)

  # Stiffness matrix
  # Material component
  Kₘ = B' * ∂𝕊∂𝔼 * B * vol
  
  # Geometric component
  auxiliar_goe_matrix = funder' * S * funder  * vol 
  Kᵧ = zeros(12,12) 
  for i=1:4
    for j=1:4
      Kᵧ[(i-1)*3+1 , (j-1) * 3 + 1] = auxiliar_goe_matrix[i,j]
      Kᵧ[(i-1)*3+2 , (j-1) * 3 + 2] = auxiliar_goe_matrix[i,j]
      Kᵧ[(i-1)*3+3 , (j-1) * 3 + 3] = auxiliar_goe_matrix[i,j]
    end
  end

  Kₑ = Kₘ + Kᵧ

end

_is_symetric(A) = (all(isapprox.(A - A', 0; rtol=1e-10)))

function matrix2vogit(𝕋::AbstractMatrix, α::Real=1)

  _is_symetric(A) || throw(ArgumentError("Tensor is not symetric"))

   v = [𝕋(1,1),𝕋(2,2),𝕋(3,3),α*𝕋(2,3), α*𝕋(1,3), α*𝕋(1,2)]' 
end

function computes_cosserat_tensor(m::SVK)

  λ, G = lambda(m)

  𝕊 = λ * tr(𝔼) * eye(3) + 2 * G * 𝔼

  ∂𝕊∂𝔼 = zeros(6, 6)
  ∂𝕊∂𝔼[1:3, 1:3] = λ * ones(3) + 2 * G * eye(3)
  ∂𝕊∂𝔼[4:6, 4:6] = G * eye(3)
  return 𝕊,∂𝕊∂𝔼 
end


function compute_B_mat(deriv::AbstractMatrix , F::AbstractMatrix)

  B = zeros(6, 12) 

  B[1:3, :] = [diagm(deriv[:,1])*F' diagm(deriv[:,2])*F' diagm(deriv[:,3])*F' diagm(deriv[:,4])*F']

  for k in 1:4
      B[4:6 , (k-1)*3 .+ (1:3)] = [ deriv[2,k] * F[:,3]' + deriv[3,k]*F[:,2]'
                                    deriv[1,k] * F[:,3]' + deriv[3,k]*F[:,1]'
                                    deriv[1,k] * F[:,2]' + deriv[2,k]*F[:,1]' ] 
  end
  return B
end
