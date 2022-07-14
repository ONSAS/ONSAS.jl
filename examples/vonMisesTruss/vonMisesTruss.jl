## Von Mises truss example problem

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

steel1  = Material( "LinearElastic", [ E1 ] )
steel2  = Material( "LinearElastic", [ E2 ] )

seccion_barra = CrossSection( "square", [sqrt(A)] )

elementGeom = ElementGeometry( "truss", seccion_barra)

# The coordinate matrix is given by
nodalCoords = [  0.  0.  0. ;
                 d   0.  h  ;
                2d   0.  0. ]
