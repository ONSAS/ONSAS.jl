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
boundary_conditions = Dict(
    "fixed" => bc₁,
    "load" => bc₂
)
# -------------------------------
# Create mesh
# -------------------------------
## Nodes
n₁ = Node((0.0, 0.0, 0.0))
n₂ = Node((d, h, 0.0))
n₃ = Node((2d, 0.0, 0.0))
nodes = [n₁, n₂, n₃]
## Connectivity
conec_nodes = Dict{String,Set}(
    "fixed" => Set([n₁, n₂]),
    "load" => Set([n₃])
)
conec_elems = Dict{String,Set}(
    "truss_s₁" => Set([[n₁, n₂]]),
    "truss_s₂" => Set([[n₂, n₃]])
)


mesh = Mesh(nodes, conec_nodes, conec_elems)
