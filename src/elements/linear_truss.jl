"""
function to linear truss element
"""

function linear_truss(material, geometry, nodalCoords::Matrix, u::Vector)

  E = material.constitutive_params[1]
  A = geometry.area

  # Element geometry
  diff = SVector{3}(nodalCoords[2, :]) - SVector{3}(nodalCoords[1, :])
  length = sqrt(diff' * diff) # Element length
  (c, s) = (diff[1], diff[3]) ./ length
  # Rotation matrix
  Qloc2glo = @SMatrix [c -s 0 0
    s c 0 0
    0 0 c -s
    0 0 s c]
  # Stiffness matrix in local system
  Kloc = E * A / length * @SMatrix [1 0 -1 0
    0 0 0 0
    -1 0 1 0
    0 0 0 0]
  # Stiffness matrix in global system
  Kglo = Qloc2glo * Kloc * Qloc2glo'
  # Internal forces
  fᵢₙₜ = Kglo * u
  # Store
  force_vectors = [fᵢₙₜ]
  tangent_matrices = [Kglo]

  return force_vectors, tangent_matrices
end