"""
function to linear truss element
"""
function linear_truss( material, geometry, nodalCoords, u  )

  E = material.constitutive_params[1]
  A = geometry.area

  diff = nodalCoords[ 2, : ] - nodalCoords[ 1, : ]

  length   = sqrt( diff' * diff )
  c        = diff[1] / length 
  s        = diff[3] / length 

  Qloc2glo = [  c -s 0  0 ;
                s  c 0  0 ;
                0 0  c -s ;
                0 0  s  c ] ;

  Kloc     = E*A/length * [ 1  0 -1 0 ;
                            0  0  0 0 ;
                            -1 0  1 0 ;
                            0  0  0 0 ] ;

  Kglo     = Qloc2glo * Kloc * transpose(Qloc2glo)

  fint = Kglo * u

  
  force_vectors    = [ fint ]
  tangent_matrices = [ Kglo ]

  return force_vectors, tangent_matrices

end