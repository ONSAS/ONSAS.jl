# Include `create_mesh` function.
include(joinpath(pkgdir(ONSAS), "examples", "cylinder_internal_pressure", "cylinder_mesh.jl"))

"""
Cylinder with internal pressure GMSH mesh and `IsotropicLinearElastic` material.

`ms` is the refinement factor of the mesh.
"""
function linear_cylinder_structure(; ms::Real=0.5)

    ## scalar parameters (dimensions in mm an MPa)
    Lₖ = 30.0       # cylinder length in 𝐞ₖ mm
    Rᵢ = 100.0      # inner radius in mm
    Rₑ = 200.0      # outer radius in mm

    E = 210.0       # Young modulus in MPa
    ν = 0.3         # Poisson ratio

    p = 25.0        # pressure in MPa
    pressure(t) = -p * t

    # -------------------------------
    # Physical entities labels
    # -------------------------------
    # material
    mat_label = "mat"
    # entities
    node_label = "node"
    faces_label = "triangle"
    elements_label = "tetrahedron"
    entities_labels = [node_label, faces_label, elements_label]
    # boundary conditions
    bc₁_label = "fixed-ui"
    bc₂_label = "fixed-uj"
    bc₃_label = "fixed-uk"
    bc₄_label = "pressure"
    bc_labels = [bc₁_label, bc₂_label, bc₃_label, bc₄_label]
    # mesh labels
    labels = [mat_label, entities_labels, bc_labels]
    # -------------------------------
    # Entities
    # -------------------------------
    # Entities types without assigned nodes, faces and elements
    vfaces = [TriangularFace(faces_label)]
    velems = [Tetrahedron(elements_label)]
    entities = StructuralEntity(velems, vfaces)
    # -------------------------------
    # Mesh
    # -------------------------------
    filename = basename(tempname())
    labels = [mat_label, entities_labels, bc_labels]
    dir = joinpath(pkgdir(ONSAS), "benchmark", "linear_cylinder_internal_pressure")
    filename = basename(tempname())
    local msh_path
    out = @capture_out begin
        msh_path = create_cylinder_mesh(Rᵢ, Rₑ, Lₖ, labels, filename, ms, dir)
    end
    gmsh_println(out)
    msh_mesh = MshFile(msh_path)
    mesh = Mesh(msh_mesh, entities)
    # Dofs
    #--------------------------------
    dof_dim = 3
    dof_u_symbol = :u
    set_dofs!(mesh, dof_u_symbol, dof_dim)
    # -------------------------------
    # Boundary conditions
    # -------------------------------
    # Dirichlet boundary conditions
    bc₁ = FixedField(:u, [1], bc₁_label)
    bc₂ = FixedField(:u, [2], bc₂_label)
    bc₃ = FixedField(:u, [3], bc₃_label)
    # Neumann boundary conditions
    bc₄ = Pressure(:u, pressure, bc₄_label)
    boundary_conditions = StructuralBoundaryCondition(bc₁, bc₂, bc₃, bc₄)
    # Assign boundary conditions to the ones defined in the mesh
    apply!(boundary_conditions, mesh)
    # -------------------------------
    # Materials
    # -------------------------------
    material = IsotropicLinearElastic(E, ν, mat_label)
    materials = StructuralMaterial(material)
    apply!(materials, mesh)
    # -------------------------------
    # Structure
    # -------------------------------

    # Create mesh and retrieve the Structure
    Structure(mesh, materials, boundary_conditions)
end;

"""
Hyper rectangle starting at `O` pint and + [Lᵢ,Lⱼ,Lₖ] to evaluate the solution with `NPOINTS` in each axis.
"""
function point_eval_handler(structure::Structure;
                            NPOINTS::Int=10)

    ## scalar parameters (dimensions in mm an MPa)
    Lₖ = 30.0                         # cylinder length in 𝐞ₖ mm
    Rₑ = 200.0                       # outer radius in mm
    Lᵢ = Lⱼ = 2.25Rₑ                # hyper rectangle origin in 𝐞ᵢ,𝐞ⱼ and  𝐞ₖ in mm
    O = (x=-Lᵢ / 2, y=-Lⱼ / 2, z=0.0)    # hyper rectangle origin in 𝐞ᵢ,𝐞ⱼ and  𝐞ₖ in mm

    # Create an hyper rectangle Lᵢ x Lⱼ x Lₖ
    x = LinRange(O.x, O.x + Lᵢ, NPOINTS)
    y = LinRange(O.y, O.y + Lⱼ, NPOINTS)
    z = LinRange(O.z, O.z + Lₖ, NPOINTS)

    vec_points = [Point(p...) for p in vec(collect(Iterators.product(x, y, z)))]
    PointEvalHandler(mesh(structure), vec_points)
end;
