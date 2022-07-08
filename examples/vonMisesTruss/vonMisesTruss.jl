
using ONSAS

## scalar parameters
E  = 2e11  # Young modulus in Pa
A  = 5e-3  # Cross-section area in m^2
d  = 1.0   # horizontal distance in m
h  = 1.0   # height in m
Fx = 0     # horizontal load in N
Fy = -1e3  # vertical   load in N

## set structs
E1 = E
E2 = E

steel1  = MaterialsData( E1, 0.3 )
steel2  = MaterialsData( E2, 0.3 )
section = ElementsData("truss", A)

# The coordinate matrix is given by
nodalCoords = [  0.  0.  0. ;
                 d   0.  h  ;
                2d   0.  0. ]

# matrix with MEBI indexes of each element (on each row)
# MEBI = Materials / Elements / BoundaryConditions / InitialConditions
MEBIValsMat = [ 0 1 1 0  ; # no material / first element / first BC / no IC
0 1 2 0  ;
1 2 0 0  ;
2 2 0 0 ] ;

# the connectivity matrix is given by a vector of vectors.
elemNodalConnec = [ [ 1 ],  [ 2 ], [ 3 ], [ 1, 2],  [ 2, 3] ]
MEBIVec         = [     1,      2,     1,       3,       4  ] 

BCs = [ BoundaryCondsData([ 1, 3, 5],[0, 0, 0],[     ], [ ]       ) ,
        BoundaryCondsData([ 3    ]  ,[0      ],[1, 5 ], [Fx, Fy ] ) ]

# assemble the sitiffness equation
KG, FG, Kmatrices, neumDofs = assembler( nodalCoords, elemNodalConnec, MEBIValsMat, MEBIVec, [steel1, steel2], [[],section], BCs )

print("KG:\n")
display(KG)
display(FG)
print("neumDofs", neumDofs)
print("\n")

KGred = copy( KG )
FGred = copy( FG )

KGred = KGred[ neumDofs, neumDofs ]
FGred = FGred[ neumDofs           ]

print("KGred", KGred,"\n")
# the system is solved.
UGred = KGred \ FGred 

UG = zeros( size( FG ) )
UG[  neumDofs ] = UGred 
