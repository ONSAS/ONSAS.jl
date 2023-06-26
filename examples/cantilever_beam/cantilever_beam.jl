using ONSAS

# parametros
b = 0.2 # width
h = 0.3 # heigth
pp = 2400

# carga distribuida
q = b * h * pp * 9.81

# -----------------------------------------
# importacion de malla
L = 3
num_elems = 4
x_coords = range(0, L, num_elems + 1)
# -----------------------------------------

# Mesh

nodes = [Node(xi, 0.0, 0.0) for xi in x_coords]

S = Rectangle(b, h) # cm^2
frames = [Frame(nodes[j], nodes[j + 1], S) for j in 1:(length(nodes) - 1)]

msh = Mesh(; nodes=nodes, elements=frames)
# -----------------------------------------

set_dofs!(msh, :u, 3) # rename create_dof
set_dofs!(msh, :θ, 3) # rename create_dof

# Materials
i = IsotropicLinearElastic(210e9, 0.3)
mat = StructuralMaterial(i => frames)

# Boundary conditions
bc1 = FixedDof(:u, [1, 2, 3])
bc2 = FixedDof(:θ, [1])

# bc3 = GlobalLoad(:u, t -> [0, 0, Py])
# bc4 = GlobalLoad(:θ, t -> [0, -Pz, 0])

bc = StructuralBoundaryCondition(bc1 => [nodes[1], nodes[end]], bc2 => [nodes[1], nodes[end]])

# Structure
s = Structure(msh, mat, bc)

# Analysis
anali = LinearStaticAnalysis(s; NSTEPS=10)

solve!(anali)
