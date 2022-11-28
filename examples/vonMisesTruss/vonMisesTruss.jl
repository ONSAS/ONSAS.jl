## Von Mises truss example problem

using ONSAS

## scalar parameters
E  = 2e11  # Young modulus in Pa
A  = 5e-3  # Cross-section area in m^2
d  = 1.0   # horizontal distance in m
h  = 1.0   # height in m
# Fx = 0     # horizontal load in N
# Fy = -1e3  # vertical   load in N

## set structs
E1 = E
E2 = E

# -------------------------------
# Materials
steel1  = Material( "LinearElastic", [ E1 ] )
steel2  = Material( "LinearElastic", [ E2 ] )

materials = [steel1, steel2]
# -------------------------------

# -------------------------------
# Geometries
node_geometry  = Geometry( "node" )
square_section = CrossSection( "square", width_y = sqrt(A) )
truss_geometry = Geometry( "truss", square_section )

geometries = [ node_geometry, truss_geometry ]
# -------------------------------

# -------------------------------
# BoundaryConditions
fixed_support    = BoundaryCondition( [1,3,5], zeros(3), "" )
load_and_support = BoundaryCondition( [3    ], [0.0   ], "my_nodal_load" )

function my_nodal_load( solution::ModelSolution, properties::ModelProperties)
    num_nodes = size( properties.mesh.nodes_coords, 1 )
    f_ext = zeros( 6*num_nodes )
    f_ext[6+5] = - 1e3 * solution.time
    return f_ext
end

boundary_conditions = [ fixed_support, load_and_support ]
# -------------------------------


# -------------------------------
# InitialConditions
initial_conditions = []
# -------------------------------


# -------------------------------
# Mesh

# The coordinate matrix is given by
nodal_coords = [  0.  0.  0. ;
                 d   0.  h  ;
                2d   0.  0. ]

elem_nodal_connec = [ [ 1 ],  [ 2 ], [ 3 ], [ 1, 2],  [ 2, 3] ]

# matrix with MGBI indexes of each element (on each row)
MGBIValsMat = [ 0 1 1 0  ; # no material / first geometry / first BC / no IC
                0 1 2 0  ;
                1 2 0 0  ;
                2 2 0 0 ]

MGBIVec     = [ 1, 2, 1, 3, 4 ]

my_mesh = Mesh( nodal_coords, elem_nodal_connec, MGBIValsMat, MGBIVec )
# -------------------------------



#fixed_node = Node( fixed_support )
#fixed_node = Node( fixed_support )
#loaded_node = Truss( steel2, square_section, load_and_support )

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
analysis_settings = AnalysisSettings( "newton_raphson", 1.0, 2.0 );
# -------------------------------

#print("type truss:", typeof(Truss2)==Symbol("Truss") )

#mallab = lector_msh( materiales, cross_section, boundary_conditions, initial_conditions )

initial_solution, model_properties = ONSAS_init( materials, geometries, boundary_conditions, initial_conditions, my_mesh, analysis_settings )

print("initial solution", initial_solution)

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
