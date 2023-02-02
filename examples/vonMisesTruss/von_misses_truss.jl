## Von Mises truss example problem
using ONSAS
using StaticArrays: SVector, SMatrix
using SparseArrays
## scalar parameters
E = 2e11  # Young modulus in Pa
ν = 0.0  # Poisson's modulus
A = 5e-3  # Cross-section area in m^2
ang = 65 # truss angle in degrees
L = 2 # Length in m 
d = L * cos(deg2rad(65))   # vertical distance in m
h = L * sin(deg2rad(65))
# Fx = 0     # horizontal load in N
Fⱼ = -3e8  # vertical   load in N
# -------------------------------
# Materials
# -------------------------------
steel = SVK(E, ν)
aluminum = SVK(E / 3, ν)
materials_dict = Dict("steel" => steel, "aluminum" => aluminum)
materials = StructuralMaterials(materials_dict)
# -------------------------------
# Geometries
# -------------------------------
## Cross section
a = sqrt(A)
s₁ = Square(a)
s₂ = Square(2a)
# -------------------------------
# Elements
# -------------------------------
elements_dict = Dict{String,AbstractElement}(
    "truss_s₁" => Truss(s₁),
    "truss_s₂" => Truss(s₂)
)
elements = StructuralElements(elements_dict)
# -------------------------------
# Boundary conditions
# -------------------------------
bc₁ = FixedDisplacementBoundaryCondition()
bc₂ = FⱼLoadBoundaryCondition(Fⱼ)
bcs_dict = Dict(
    "fixed" => bc₁,
    "load" => bc₂
)
boundary_conditions = StructuralBoundaryConditions(bcs_dict)
# -------------------------------
# Create mesh
# -------------------------------
## Nodes
n₁ = Node((0.0, 0.0, 0.0))
n₂ = Node((d, h, 0.0))
n₃ = Node((2d, 0.0, 0.0))
vec_nodes = [n₁, n₂, n₃]
## Elements connectivity
elem₁_nodes = [n₁, n₂]
elem₂_nodes = [n₂, n₃]
vec_conec_elems = [elem₁_nodes, elem₂_nodes]
mesh = Mesh(vec_nodes, vec_conec_elems)

# -------------------------------
# Apply MEBI
# -------------------------------
# Materials
material_sets = Dict{String,Set{Int}}(
    "steel" => Set{Int}([1]),
    "aluminum" => Set{Int}([2]),
)
add_set!(materials, material_sets)
# Elements
element_sets = Dict{String,Set{Int}}(
    "truss_s₁" => Set{Int}([1]),
    "truss_s₂" => Set{Int}([2]),
)
add_set!(elements, element_sets)
# Boundary conditions
bc_sets = Dict{String,Set{Int}}(
    "fixed" => Set{Int}([1, 3]),
    "load" => Set{Int}([2])
)
add_set!(boundary_conditions, bc_sets)
# -------------------------------
# Create Structure
# -------------------------------
s = Structure(mesh, materials, elements, boundary_conditions)


