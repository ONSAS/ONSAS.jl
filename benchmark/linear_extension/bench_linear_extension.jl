# Include `create_mesh` function.
include(joinpath(pkgdir(ONSAS), "examples", "linear_extension", "linear_extension_mesh.jl"))

"""
Linear extension  GMSH mesh and `IsotropicLinearElastic` material.

`ms` is the refinement factor of the mesh.
"""
function linear_extension_structure(; ms=0.5)

    # x, y and z dimensions of the box in the mesh respectively.
    Lᵢ = 2.0                         # Dimension in x of the box in m 
    Lⱼ = 1.0                         # Dimension in y of the box in m
    Lₖ = 1.0                         # Dimension in z of the box in m

    ## scalar parameters
    E = 2.0                          # Young modulus in Pa
    ν = 0.4                          # Poisson's ratio
    mat_label = "mat"
    mat = IsotropicLinearElastic(E, ν, mat_label)

    # Tension load in Pa.
    p = 3

    # Material types without assigned elements.
    materials = StructuralMaterials(mat)

    # Dirichlet boundary conditions 
    bc₁_label = "fixed-ux"
    bc₂_label = "fixed-uj"
    bc₃_label = "fixed-uk"
    bc₄_label = "tension"
    bc₁ = FixedDofBoundaryCondition([:u], [1], bc₁_label)
    bc₂ = FixedDofBoundaryCondition([:u], [2], bc₂_label)
    bc₃ = FixedDofBoundaryCondition([:u], [3], bc₃_label)

    # Neumann boundary conditions 
    bc₄ = GlobalLoadBoundaryCondition([:u], t -> [tension(t), 0, 0], bc₄_label)

    boundary_conditions = StructuralBoundaryConditions(bc₁, bc₂, bc₃, bc₄)
    bc_labels = [bc₁_label, bc₂_label, bc₃_label, bc₄_label]

    # Entities types without assigned nodes, faces and elements
    faces_label = "triangle"
    elems_label = "tetrahedron"
    vfaces = [TriangularFace(faces_label)]
    velems = [Tetrahedron(elems_label)]
    entities_labels = [faces_label, elems_label]
    entities = StructuralEntities(velems, vfaces)

    # Create mesh and retrieve the Structure

    filename = basename(tempname())
    labels = [mat_label, entities_labels, bc_labels]
    dir = joinpath(pkgdir(ONSAS), "benchmark", "linear_extension")
    msh_file = MshFile(create_linear_extension_mesh(Lᵢ, Lⱼ, Lₖ, labels, filename, ms, dir))

    Structure(msh_file, materials, boundary_conditions, entities)
end
