# --------------------------------------------------
# Cylinder submitted to an Internal Pressure Example  
#----------------------------------------------------
using ONSAS.StaticAnalyses
using Test: @test, @testset
## scalar parameters (dimensions in mm an MPa)
const p = 1e-2 # internal pressure in MPa
const L = 0.75 # cylinder length in 𝐞ₖ mm
const Rᵢ = 100 # inner radius in mm
const Rₑ = 200 # outer radius in mm
const E = 210  # Young modulus in MPa
const ν = 0.3  # Poisson ratio
# ------------------------------------------
# Case 1 - `IsoptrpicLinearElastic` material
# -------------------------------------------
# -------------------------------
# Entities
# -------------------------------
# Entities types without assigned nodes, faces and elements
faces_label = "triangle"
elements_label = "tetrahedron"
vfaces = [TriangularFace(faces_label)]
velems = [Tetrahedron(elements_label)]
s_entities = StructuralEntities(velems, vfaces)
entities_labels = [faces_label, elements_label]
# -------------------------------
# Mesh
# -------------------------------
filename = "cylinder"
msh_file = MshFile(file_name_mesh)
s_mesh = Mesh(mesh_file, s_entities)


#=
# -------------------------------
# Boundary conditions
# -------------------------------
# Fixed dofs
bc₁_label = "fixed-ux"
bc₁ = FixedDofBoundaryCondition([:u], [1], bc₁_label)
bc₂_label = "fixed-uj"
bc₂ = FixedDofBoundaryCondition([:u], [2], bc₂_label)
bc₃_label = "fixed-uk"
bc₃ = FixedDofBoundaryCondition([:u], [3], bc₃_label)
# Load
bc₄_label = "pressure"
bc₄_pressure_function(t) = -p * t # pressure function towards -𝐞ᵣ
bc₄ = LocalPressureBoundaryCondition([:u], t -> [bc₄_pressure_function(t)], bc₄_label)
# Create the  bcs vector
bc_labels = [bc₁_label, bc₂_label, bc₃_label, bc₄_label]
# BoundaryConditions types without assigned node, feces and elements
bcs = [bc₁, bc₂, bc₃, bc₄]
s_boundary_conditions = StructuralBoundaryConditions(bcs)
# Assign boundary conditions to the ones defined in the mesh
apply!(s_boundary_conditions, s_mesh)
# -------------------------------
# Materials
# -------------------------------
mat_label = "linear_elastic"
mat = IsotropicLinearElastic(E, ν, mat_label)
s_materials = StructuralMaterials([mat])
# Assign material to the elements defined in the mesh
apply!(s_materials, s_mesh)
=#