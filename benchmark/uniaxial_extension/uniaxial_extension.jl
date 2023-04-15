using LinearAlgebra

# Strain energy function of an hyperelastic material used in the benchmark.
strain_energy_svk(𝔼::AbstractMatrix, λ::Real, G::Real) = (λ / 2) * tr(𝔼)^2 + G * tr(𝔼^2)

# Include `create_mesh` function.
include(joinpath(pkgdir(ONSAS), "examples", "uniaxial_extension", "uniaxial_cube_mesh.jl"))

"""
Uniaxial extension Case 2 - GMSH mesh and `HyperElastic` material.

`ms` is the refinement factor of the mesh.
"""
function uniaxial_extension_structure(; ms=0.5)
    # x, y and z dimensions of the box in the mesh respectively.
    Lᵢ = 2.0
    Lⱼ = 1.0
    Lₖ = 1.0
    # Young's modulus in Pa.
    E = 1.0
    # Poisson's ratio.
    ν = 0.3
    mat_label = "svkHyper"
    svk = SVK(E, ν, mat_label)
    λ, G = lame_parameters(svk)
    svk_hyper_elastic = HyperElastic([λ, G], strain_energy_svk, "svkHyper")
    # Tension load in Pa.
    p = 3

    # Material types without assigned elements.
    materials = StructuralMaterials(svk_hyper_elastic)

    # Redefine the load boundary condition.
    bc₁_label = "fixed-ux"
    bc₂_label = "fixed-uj"
    bc₃_label = "fixed-uk"
    bc₄_label = "tension"
    bc₁ = FixedDofBoundaryCondition([:u], [1], bc₁_label)
    bc₂ = FixedDofBoundaryCondition([:u], [2], bc₂_label)
    bc₃ = FixedDofBoundaryCondition([:u], [3], bc₃_label)
    bc₄ = LocalPressureBoundaryCondition([:u], t -> [p * t], bc₄_label)
    bc_labels = [bc₁_label, bc₂_label, bc₃_label, bc₄_label]

    # BoundaryConditions types without assigned node, feces and elements.
    boundary_conditions = StructuralBoundaryConditions(bc₁, bc₂, bc₃, bc₄)

    # Entities types without assigned nodes, faces and elements.
    faces_label = "triangle"
    elems_label = "tetrahedron"
    vfaces = [TriangularFace(faces_label)]
    velems = [Tetrahedron(elems_label)]
    entities = StructuralEntities(velems, vfaces)

    entities_labels = [faces_label, elems_label]
    filename = basename(tempname())
    labels = [mat_label, entities_labels, bc_labels]
    mesh = MshFile(create_mesh(Lᵢ, Lⱼ, Lₖ, labels, filename, ms))
    Structure(mesh, materials, boundary_conditions, entities)
end
