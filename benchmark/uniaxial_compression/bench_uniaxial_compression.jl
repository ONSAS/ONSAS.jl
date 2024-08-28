using LinearAlgebra

# Strain energy function of an hyperelastic material used in the benchmark.
function strain_energy_neo(𝔼::AbstractMatrix, K::Real, μ::Real)
    # Right hand Cauchy strain tensor
    ℂ = Symmetric(2 * 𝔼 + eye(3))
    J = sqrt(det(ℂ))
    # First invariant
    I₁ = tr(ℂ)
    # Strain energy function
    Ψ = μ / 2 * (I₁ - 2 * log(J)) + K / 2 * (J - 1)^2
end

# Include `create_mesh` function.
include(joinpath(pkgdir(ONSAS), "examples", "uniaxial_extension", "uniaxial_mesh.jl"))

"""
Uniaxial compression Case 2 - GMSH mesh and `HyperElastic` material.

`ms` is the refinement factor of the mesh.
"""
function uniaxial_compression_structure(; ms=0.5)

    # x, y and z dimensions of the box in the mesh respectively.
    Lᵢ = 2.0
    Lⱼ = 1.0
    Lₖ = 1.0

    E = 1.0                    # Young modulus in Pa
    ν = 0.3                    # Poisson's ratio
    K = E / (3 * (1 - 2 * ν))  # Bulk modulus in Pa
    μ = G = E / (2 * (1 + ν))  # Second Lamé parameter in Pa
    mat_label = "neoHyper"
    neo_hookean_hyper = HyperElastic([K, μ], strain_energy_neo, "neoHyper")

    # Tension load in Pa.
    p = -1

    # Material types without assigned elements.
    materials = StructuralMaterial(neo_hookean_hyper)

    # Dirichlet boundary conditions
    bc₁_label = "fixed-ux"
    bc₂_label = "fixed-uj"
    bc₃_label = "fixed-uk"
    bc₄_label = "tension"
    bc₁ = FixedField(:u, [1], bc₁_label)
    bc₂ = FixedField(:u, [2], bc₂_label)
    bc₃ = FixedField(:u, [3], bc₃_label)

    # Neumann boundary conditions
    bc₄ = Pressure(:u, t -> p * t, bc₄_label)
    bc_labels = [bc₁_label, bc₂_label, bc₃_label, bc₄_label]

    # BoundaryConditions types without assigned node, feces and elements.
    boundary_conditions = StructuralBoundaryCondition(bc₁, bc₂, bc₃, bc₄)

    # Entities types without assigned nodes, faces and elements.
    faces_label = "triangle"
    elems_label = "tetrahedron"
    vfaces = [TriangularFace(faces_label)]
    velems = [Tetrahedron(elems_label)]
    entities = StructuralEntity(velems, vfaces)
    entities_labels = [faces_label, elems_label]

    # Create mesh and retrieve the Structure
    entities_labels = [faces_label, elems_label]
    filename = basename(tempname())
    labels = [mat_label, entities_labels, bc_labels]
    dir = joinpath(pkgdir(ONSAS), "benchmark", "uniaxial_compression")
    mesh = MshFile(create_uniaxial_mesh(Lᵢ, Lⱼ, Lₖ, labels, filename, ms, dir))
    return Structure(mesh, materials, boundary_conditions, entities)
end
