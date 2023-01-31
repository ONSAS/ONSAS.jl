## Von Mises truss example problem

using ONSAS

## scalar parameters
E = 2e11  # Young modulus in Pa
ν = 0.0  # Poisson's modulus
A = 5e-3  # Cross-section area in m^2
ang = 65 # truss angle in degrees
L = 2 # Length in m 
d = L * cos(deg2rad(65))   # vertical distance in m
h = L * sin(deg2rad(65))
d = 0.0  # horizontal distance in m
# Fx = 0     # horizontal load in N
Fⱼ = -3e8  # vertical   load in N
# -------------------------------
# Materials
# -------------------------------
steel = SVK(E, ν)
materials = [steel]
# -------------------------------
# Geometries
# -------------------------------
## Cross section
dim = sqrt(A)
s = Square(dim)
## Nodes
n₁ = Node((0.0, 0.0))
n₂ = Node((d, h))
n₃ = Node((2d, 0.0))
# -------------------------------
# Boundary conditions
# -------------------------------
bc₁ = FixedDisplacementBoundaryCondition()
bc₂ = FⱼLoadBoundaryCondition(Fⱼ)
push!([n₁, n₃], bc₁)
push!([n₂], bc₂)
# -------------------------------
# Elements
# -------------------------------
t₁ = Truss([n₁, n₂], steel, s)
# -------------------------------
# Test internal force
# -------------------------------
u_e = [[0, 0], [0, 0]]
@show K₁ = stiffness_matrix(t₁, u_e)
@show f₁ = internal_force(t₁, u_e)



#=


Section = Rectangle(dim, dim)
Square_section = CrossSection(Section)
## Geometries
node_geometry = Geometry(Node())
truss_geometry = Geometry(Truss(), Square_section)

Geometries = [node_geometry, truss_geometry]
# -------------------------------

# -------------------------------
# BoundaryConditions
fixed_support = DispsBoundaryCondition([1, 3, 5], zeros(3))
y_support = DispsBoundaryCondition([3], [0.0])
# load_and_support = DispsBoundaryCondition([3], [0.0], "my_nodal_load")

function my_nodal_load(solution::ModelSolution, properties::ModelProperties)
    num_nodes = size(properties.mesh.nodes_coords, 1)
    f_ext = zeros(6 * num_nodes)
    f_ext[6+5] = -1e3 * solution.time
    return f_ext
end

Fz = -1e3
loadsBaseVals = [0, 0, 0, 0, -1e3, 0]
loadsCoordSystem = "Global"
LoadsBC = [LoadsBoundaryCondition(loadsBaseVals, loadsCoordSystem)]

DofsBC = [fixed_support, y_support]
# -------------------------------


# -------------------------------
# InitialConditions
IC = []
# -------------------------------


# -------------------------------
# Mesh

# The coordinate matrix is given by
nodal_coords = [0.0 0.0 0.0
    d 0.0 h
    2d 0.0 0.0]

elem_nodal_connec = [[1], [2], [3], [1, 2], [2, 3]]

# matrix with MGBI indexes of each element (on each row)
MGBIValsMat = [0 1 0 1 0 # no material / first geometry / load BC / supportBC /no IC
    0 1 1 2 0
    1 2 0 0 0
    2 2 0 0 0]

MGBIVec = [1, 2, 1, 3, 4]

StrMesh = Mesh(nodal_coords, elem_nodal_connec, MGBIValsMat, MGBIVec)
# -------------------------------
StrConvergenceSettings = ConvergenceSettings()
Algorithm = NewtonRaphson(1.0, 1.0)

initial_solution, model_properties = ONSAS_init(Materials, Geometries, LoadsBC, DofsBC, IC, StrMesh, StrConvergenceSettings, Algorithm)
stop
print("initial solution \n", initial_solution)


#fixed_node = Node( fixed_support )
#fixed_node = Node( fixed_support )
#loaded_node = Truss( steel2, Square_section, load_and_support )

# elementos = AbstractElement[]

# indices = [1,2]

# for j in indices
#     print("j ", j,"\n")
#     push!(elementos, Node( fixed_support ) )
#     print(elementos,"\n")
# end

# print("\n\nvalor ",  elementos[1],"\n")
# print("\n\nTYPEOF ", typeof( elementos[1]),"\n")

# elementos[1].connectivity = [ 4]
# print("PRUEBA ",  elementos,"\n")

# elementos[2].connectivity = [ 3]
# print("PRUEBA 2 ",  elementos,"\n")

# mallab = MeshB( nodal_coords, elementos)

# print( "\n\n mallabb", mallab,"\n")
#print(" inertia: ", mallab.elements[2].cross_section.inertia_y )

# -------------------------------
# AnalysisSettings
# analysis_settings = AnalysisSettings("newton_raphson", 1.0, 2.0);
Analysis_settings = AnalysisSettings()

# -------------------------------

#print("type truss:", typeof(Truss2)==Symbol("Truss") )

#mallab = lector_msh( materiales, cross_section, boundary_conditions, initial_conditions )

initial_solution, model_properties = ONSAS_init(materials, geometries, boundary_conditions, initial_conditions, my_mesh, Analysis_settings)

print("initial solution \n", initial_solution)

#solutions = ONSAS_solve( model_properties, initial_solution, verbosity=true )


# print("KG:\n")
# display(KG)
# display(FG)
# print("neumDofs", neumDofs)
# print("\n")

# KGred = copy( KG )
# FGred = copy( FG )

# KGred = KGred[ neumDofs, neumDofs ]
# FGred = FGred[ neumDofs           ]

# print("KGred", KGred,"\n")
# # the system is solved.
# UGred = KGred \ FGred

# UG = zeros( size( FG ) )
# UG[  neumDofs ] = UGred

=#