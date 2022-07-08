# # Uniaxial extension problem
#

using FEMAssembler

# scalar parameters
E  = 1.0 # Young modulus in Pa
ν  = .3  # Poisson's ratio
p  = 2   # Tension load in Pa
Lx = 2   # Dimension in x of the boxD in m 
Ly = 1   # Dimension in y of the boxD in m
Lz = 1   # Dimension in z of the boxD in m
 
# Lamé parameters
λ = E*ν / ( (1+ν) * (1-2*ν) ) 
μ = E / ( 2 * (1+ν) ) 

# set material structs
SVKmat = MaterialsData( E, ν );

# set element structs
triangl = ElementsData( "face", 0 )
tetra   = ElementsData( "tetrahedron" , 0 )

# Set the loads
loadsTimeFact= t-> p*t 

# boundary conditions
BCs = [ BoundaryCondsData( [ ], [ ], [1], [ loadsTimeFact( 0 ) ] ),  # tension
        BoundaryCondsData( [1], [0], [ ], [ ] ),  # non-friction constraint x=0
        BoundaryCondsData( [3], [0], [ ], [ ] ),  # non-friction constraint y=0
        BoundaryCondsData( [5], [0], [ ], [ ] ) ] # non-friction constraint z=0

# the connectivity matrix is given by a vector of vectors.
elemNodalConnec =[ [ 5 8 6   ] ,  # loaded face
                   [ 6 8 7   ] ,  # loaded face
                   [ 4 1 2   ] ,  # x=0 supp face
                   [ 4 2 3   ] ,  # x=0 supp face
                   [ 6 2 1   ] ,  # y=0 supp face
                   [ 6 1 5   ] ,  # y=0 supp face
                   [ 1 4 5   ] ,  # z=0 supp face
                   [ 4 8 5   ] ,  # z=0 supp face
                   [ 1 4 2 6 ] ,  # tetrahedron
                   [ 6 2 3 4 ] ,  # tetrahedron
                   [ 4 3 6 7 ] ,  # tetrahedron
                   [ 4 1 5 6 ] ,  # tetrahedron
                   [ 4 6 5 8 ] ,  # tetrahedron
                   [ 4 7 6 8 ]  ] # tetrahedron

# matrix with MEBI indexes of each property (on each row)
# MEBI = Materials / Elements / BoundaryConditions / InitialConditions
MEBIValsMat = [ 0 1 1 0 ;  # BC 1
                0 1 2 0 ;  # BC 2
                0 1 3 0 ;  # BC 3
                0 1 4 0 ;  # BC 4
                1 2 0 0 ]  

MEBIVec = Int.( vcat( ones(2), 2*ones(2), 3*ones(2), 4*ones(2), 5*ones(6) ) )

# The coordinate matrix is given by
nodalCoords = [ 0    0    0 ; 
                0    0   Lz ; 
                0   Ly   Lz ; 
                0   Ly    0 ; 
                Lx   0    0 ; 
                Lx   0   Lz ; 
                Lx  Ly   Lz ; 
                Lx  Ly    0 ] 

finalTime = 1;
deltaTime = .5 ; 
                
                
print("Starting time loop.\n")
currTime = 0.0

while currTime < finalTime
    global currTime = currTime + deltaTime 
    print("Solving time: ", currTime,"\n")

    BCs[1].NeumannNodalVals = [loadsTimeFact( currTime )] 
    print( "  current pressure ", BCs[1].NeumannNodalVals, "\n" )

    # Run FEMAssembler                             
    KG, FG, Kmatrices = assembler( nodalCoords, elemNodalConnec, MEBIValsMat, MEBIVec, [SVKmat], [ triangl, tetra ], BCs, [], verbosityBool=true )

end


