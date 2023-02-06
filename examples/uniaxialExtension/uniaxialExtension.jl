## Uniaxial extension example
using ONSAS
using LinearAlgebra
using Test
## scalar parameters
# scalar parameters
E = 1.0 # Young modulus in Pa
ν = 0.3  # Poisson's ratio
p = 2   # Tension load in Pa
Lᵢ = 2.0   # Dimension in x of the box in m 
Lⱼ = 1.0   # Dimension in y of the box in m
Lₖ = 1.0   # Dimension in z of the box in m
# -------------------------------
# Create mesh
# -------------------------------
n₁ = Node((0.0, 0.0, 0.0))
n₂ = Node((0.0, 0.0, Lₖ))
n₃ = Node((0.0, Lⱼ, Lₖ))
n₄ = Node((0.0, Lⱼ, 0.0))
n₅ = Node((Lᵢ, 0.0, 0.0))
n₆ = Node((Lᵢ, 0.0, Lₖ))
n₇ = Node((Lᵢ, Lⱼ, Lₖ))
n₈ = Node((Lᵢ, Lⱼ, 0.0))
vec_nodes = [n₁, n₂, n₃, n₄, n₅, n₆, n₇, n₈]
## Elements connectivity
triangle₁_nodes = [n₅, n₈, n₆]
triangle₂_nodes = [n₆, n₈, n₇]
triangle₃_nodes = [n₄, n₁, n₂]
triangle₄_nodes = [n₄, n₂, n₃]
triangle₅_nodes = [n₆, n₂, n₁]
triangle₆_nodes = [n₆, n₁, n₅]
triangle₇_nodes = [n₁, n₄, n₅]
triangle₈_nodes = [n₄, n₈, n₅]
tetra₁_nodes = [n₁, n₄, n₂, n₆]
tetra₂_nodes = [n₆, n₂, n₃, n₄]
tetra₃_nodes = [n₄, n₃, n₆, n₇]
tetra₄_nodes = [n₄, n₁, n₅, n₆]
tetra₅_nodes = [n₄, n₆, n₅, n₈]
tetra₆_nodes = [n₄, n₇, n₆, n₈]
vec_conec_elems = [
    triangle₁_nodes, triangle₂_nodes, triangle₃_nodes, triangle₄_nodes,
    triangle₅_nodes, triangle₆_nodes, triangle₇_nodes, triangle₈_nodes,
    tetra₁_nodes, tetra₂_nodes, tetra₃_nodes,
    tetra₄_nodes, tetra₅_nodes, tetra₆_nodes]
# Mesh
mesh = Mesh(vec_nodes, vec_conec_elems)
# find elements with 3 and 4 nodes and create a elementIndexVector
tetra_indexes = ElementIndex.(findall(x -> length(x) == 4, element_nodes(mesh)))
triangle_indexes = ElementIndex.(findall(x -> length(x) == 3, element_nodes(mesh)))
# -------------------------------
# Materials
# -------------------------------
svk = SVK(E, ν)
materials_dict = Dict("mySVK" => svk)
s_materials = StructuralMaterials(materials_dict)
# Sets
m_sets = Dict(
    "mySVK" => Set(tetra_indexes),
)
add_set!(s_materials, m_sets)
# -------------------------------
# Elements
# -------------------------------
elements_dict = Dict{String,AbstractElement}(
    "triangles" => Triangle(),
    "svk_tetras" => Tetrahedron()
)
s_elements = StructuralElements(elements_dict)
# Sets
e_sets = Dict(
    "triangles" => Set(triangle_indexes),
    "svk_tetras" => Set(tetra_indexes),
)
add_set!(s_elements, e_sets)


