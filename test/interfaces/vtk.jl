using Test
using WriteVTK
using ONSAS.VTK
using ONSAS.Entities
using ONSAS.TriangularFaces
using ONSAS.Tetrahedrons
using ONSAS.Trusses
using ONSAS.Frames
using ONSAS.CrossSections
using ONSAS.Circles
using ONSAS.Squares
using ONSAS.Meshes
using ONSAS.Nodes

@testset "Writing a tetrahedron mesh to a vtk file." begin
    Lx = 10.0
    Ly = 20.0
    Lz = 30.0

    # Mesh
    n1 = Node(0.0, 0.0, 0.0)
    n2 = Node(0.0, 0.0, Lz)
    n3 = Node(0.0, Ly, Lz)
    n4 = Node(0.0, Ly, 0.0)
    n5 = Node(Lx, 0.0, 0.0)
    n6 = Node(Lx, 0.0, Lz)
    n7 = Node(Lx, Ly, Lz)
    n8 = Node(Lx, Ly, 0.0)
    vec_nodes = [n1, n2, n3, n4, n5, n6, n7, n8]
    m = Mesh(; nodes = vec_nodes)

    ## Faces
    f1 = TriangularFace(n5, n8, n6)
    f2 = TriangularFace(n6, n8, n7)
    f3 = TriangularFace(n4, n1, n2)
    f4 = TriangularFace(n4, n2, n3)
    f5 = TriangularFace(n6, n2, n1)
    f6 = TriangularFace(n6, n1, n5)
    f7 = TriangularFace(n1, n4, n5)
    f8 = TriangularFace(n4, n8, n5)
    vec_faces = [f1, f2, f3, f4, f5, f6, f7, f8]
    append!(faces(m), vec_faces)
    ## Entities
    t1 = Tetrahedron(n1, n4, n2, n6)
    t2 = Tetrahedron(n6, n2, n3, n4)
    t3 = Tetrahedron(n4, n3, n6, n7)
    t4 = Tetrahedron(n4, n1, n5, n6)
    t5 = Tetrahedron(n4, n6, n5, n8)
    t6 = Tetrahedron(n4, n7, n6, n8)
    vec_elems = [t1, t2, t3, t4, t5, t6]
    append!(elements(m), vec_elems)

    u_dim = 3
    set_dofs!(m, :u, u_dim)
    temp_dim = 1
    set_dofs!(m, :T, temp_dim)

    filename = "tetrahedron_unit_test_vtk"
    vtk_mesh = VTKMeshFile(filename, m)
    fs = close(vtk_mesh.vtk)
    @test only(fs) == "$filename.vtu"
    @test vtk_mesh.vtk.Ncls == num_elements(m)
    @test vtk_mesh.vtk.Npts == num_nodes(m)

    vec_nodal_dof_data = to_vtk(collect(1:num_dofs(m, :u)))
    scalar_nodal_dof_data = to_vtk(collect(1:num_dofs(m, :T)))
    scalar_cell_data = to_vtk(rand(num_elements(m)))
    tensor_cell_data = [rand(3, 3) for _ in elements(m)]
    VTKMeshFile(filename, m) do vtx
        write_node_data(vtx, vec_nodal_dof_data, "vectorial_nodal_data";
            component_names = ["sx", "sy", "sz"])
        write_node_data(
            vtx, scalar_nodal_dof_data, "scalar_nodal_data"; component_names = ["T"])
        write_cell_data(vtx, scalar_cell_data, "scalar_cell_data"; component_names = ["σ"])
        write_cell_data(vtx, tensor_cell_data, "tensor_cell_data";
            component_names = ["σxx", "σyy", "σzz", "τyz", "τxz", "τxy", "τzy", "τzx",
                "τyx"])
    end
end
@testset "Writing a truss mesh to a vtk file." begin
    V = 1.0
    H = 1.0
    d = 0.1
    a = d
    n1 = Node(0.0, 0.0, 0.0)
    n2 = Node(V, 0.0, H)
    n3 = Node(2V, 0.0, 0.0)
    ns = [n1, n2, n3]
    s1 = Circle(d)
    s2 = Square(a)
    truss_left = Truss(n1, n2, s1)
    truss_right = Truss(n2, n3, s2)
    es = [truss_left, truss_right]
    m = Mesh(; nodes = ns, elements = es)

    u_dim = 3
    set_dofs!(m, :u, u_dim)
    temp_dim = 1
    set_dofs!(m, :T, temp_dim)

    vec_nodal_dof_data = to_vtk(collect(1:num_dofs(m, :u)))
    scalar_nodal_dof_data = to_vtk(collect(1:num_dofs(m, :T)))
    scalar_cell_data = to_vtk(rand(num_elements(m)))
    tensor_cell_data = [rand(3, 3) for _ in elements(m)]

    filename = "truss_test_vtk"
    vtk_mesh = VTKMeshFile(filename, m)
    fs = close(vtk_mesh.vtk)
    @test only(fs) == "$filename.vtu"
    @test vtk_mesh.vtk.Ncls == num_elements(m)
    @test vtk_mesh.vtk.Npts == num_nodes(m)

    VTKMeshFile(filename, m) do vtx
        write_node_data(vtx, vec_nodal_dof_data, "vectorial_nodal_data";
            component_names = ["sx", "sy", "sz"])
        write_node_data(
            vtx, scalar_nodal_dof_data, "scalar_nodal_data"; component_names = ["T"])
        write_cell_data(vtx, scalar_cell_data, "scalar_cell_data"; component_names = ["σ"])
        write_cell_data(vtx, tensor_cell_data, "tensor_cell_data";
            component_names = ["σxx", "σyy", "σzz", "τyz", "τxz", "τxy", "τzy", "τzx",
                "τyx"])
    end
end

@testset "Writing a frame mesh to a vtk file." begin
    L = 10.0
    N = 10
    d = 0.01
    s = Circle(d)
    x_coords = range(0, L, N + 1)
    ns = [Node(promote(xi, 0, 0.0)) for xi in x_coords]  # to be safe when arguments supplied to Node command are not of  same types
    frames = [Frame(ns[j], ns[j + 1], s) for j in 1:(length(ns) - 1)]
    m = Mesh(; nodes = ns, elements = frames)
    set_dofs!(m, :u, 3)
    set_dofs!(m, :θ, 3)

    filename = "frame_test_vtk"
    vtk_mesh = VTKMeshFile(filename, m)
    fs = close(vtk_mesh.vtk)
    @test only(fs) == "$filename.vtu"
    @test vtk_mesh.vtk.Ncls == num_elements(m)
    @test vtk_mesh.vtk.Npts == num_nodes(m)

    vec_linear_dof_data = to_vtk(collect(1:num_dofs(m, :u)))
    vec_angular_dof_data = to_vtk(collect(1:num_dofs(m, :θ)))
    scalar_cell_data = to_vtk(rand(num_elements(m)))
    tensor_cell_data = [rand(3, 3) for _ in elements(m)]

    VTKMeshFile(filename, m) do vtx
        write_node_data(vtx, vec_linear_dof_data, "vec_linear_dof_data";
            component_names = ["sx", "sy", "sz"])
        write_node_data(vtx, vec_angular_dof_data, "vec_angular_dof_data";
            component_names = ["tx", "ty", "tz"])
        write_cell_data(vtx, scalar_cell_data, "scalar_cell_data"; component_names = ["σ"])
        write_cell_data(vtx, tensor_cell_data, "tensor_cell_data";
            component_names = ["σxx", "σyy", "σzz", "τyz", "τxz", "τxy", "τzy", "τzx",
                "τyx"])
    end
end
