# SPDX-License-Identifier: MIT
"""
    computeFEM2GridInterpMatrix( Nodes, Conec, intGrid, numCompon, order )

 function for computation of interpolation matrix from FEM mesh to regular grid
    from displacements in tetraedra fem grid to a regular grid mesh

  Inputs:
    - Nodes: matrix with coordinates of nodes of the FEMmesh
    - Conec: vector of vectors. at i-th entry has the vector with the connectivity indexes
    - numComponentes: number of components of the magnitude to be interpolated
    - order: 0 interpolates from elements centroids
            1 interpolates from nodes
"""
function computeFEM2GridInterpMatrix( Nodes, Conec, intGrid, numCompon, order )
         
  numVoxelsX   =       intGrid.voxelNums[1]
  numVoxelsXY  = prod( intGrid.voxelNums[1:2] )
  numVoxels    = prod( intGrid.voxelNums )

  numNodes  = size( Nodes, 1)
  numElems  = length( Conec )

  typeof(Conec) != Vector{Vector{Int}} && error("type of Conec: ", typeof(Conec) )

  booleanInterFound = zeros( Bool, numVoxels )

  rowIndexes = Vector{Int}()
  colIndexes = Vector{Int}()
  nonZerVals = Vector{Float64}()

  if order == 1
    totalRows = numCompon * numVoxels
    totalCols = numCompon * numNodes
  end

  # loop in elements for see which grid nodes are in an element.
  for i in (1:numElems)

    nodesElem   = Conec[i]
    coordesElem = Nodes[ nodesElem , :]
    
    # coordenadas de bounding box de elemento
    mins = vec( minimum( coordesElem, dims=1 ) )
    maxs = vec( maximum( coordesElem, dims=1 ) )
    
    # indexes of grid within bounding box
    indsIni, indsEnd = ranges( mins, maxs, intGrid )

    if length( indsIni ) > 0

      for kg in (indsIni[3]:indsEnd[3])
        for jg in (indsIni[2]:indsEnd[2])
          for ig in (indsIni[1]:indsEnd[1])

            # number of voxel
            ind = (kg-1) * numVoxelsXY + (jg-1)*numVoxelsX + ig
            
            # if an interpolation of the node was not already obtained
            if !(booleanInterFound[ind])

              pointCand = intGrid.startVox + intGrid.voxelWidths .* ([ ig, jg, kg ] .- 1 )

              # verify if it belongs to the current element
              itIsInTetra, shapeFuncs = checkInTetra( coordesElem, pointCand )

              # if it is in the tetrahedron
              if itIsInTetra
                # then we found it
                booleanInterFound[ind] = true

                dofsGridInd = nodes2dofs( ind, numCompon )

                # if order == 0
                #    # assign interpolation of grid node ind to element i
                #    interpMatrix[ind, i] = 1
                
                if order == 1
  
                  for indShapeFunc in (1:length(shapeFuncs))

                    dofsNode = nodes2dofs( nodesElem[indShapeFunc], numCompon )
  
                    append!( rowIndexes, dofsGridInd )
                    append!( colIndexes, dofsNode    )
                    append!( nonZerVals, shapeFuncs[ indShapeFunc ]*ones(numCompon) )

                  end

                end # if order 

              end # if it is in tetra
            end # if node already analyzed
          end # loop x index
        end # loop y index
      end # loop z index
    end # if length inds > 0

  end # loop elements

  return sparse( rowIndexes, colIndexes, nonZerVals, totalRows, totalCols )
end

"""
TO DO
"""

function checkInTetra( nodesTetra, pointCand )
 
  A          = zeros(4,4)
  A[:, 1]   .= 1.0 
  A[:, 2:4] .= nodesTetra
  volTot     = det(A) / 6.0

  volsRel = zeros(4) 
 
  boolInTetra = true
  i=0
  while boolInTetra && i<4
    i = i + 1

    B        = copy(A)
    B[i,2:4] = pointCand
    volsRel[i] = ( det( B ) / 6.0 ) / volTot ;
    
    boolInTetra = boolInTetra && ( volsRel[i] >= -1.0e-10 )
  end
  
  if boolInTetra
    shapeFuncs = volsRel
  else
    shapeFuncs = []
  end

  return boolInTetra, shapeFuncs
end
"""
    interpolate_over_grid( grid, magnitudes_on_grid, nodes_to_eval )

Function used to compute intensities at nodes given by matrix nodes_to_eval using the intensity values given by magnitudes_on_grid.
"""
function interpolate_over_grid( grid, magnitudes_on_grid, nodes_to_eval )
    
    # check use of interpolation function 
    xs = grid.startVox[1]:grid.voxelWidths[1]:grid.endVox[1]
    ys = grid.startVox[2]:grid.voxelWidths[2]:grid.endVox[2]
    zs = grid.startVox[3]:grid.voxelWidths[3]:grid.endVox[3]
    
#    grid_int_interp = scale( interpolate( magnitudes_on_grid ), (xs,ys,zs))
    
    # working but not easy to compute gradients
    grid_int_interp = linear_interpolation( (xs,ys,zs), magnitudes_on_grid, extrapolation_bc = 0.0 )

    if isa(nodes_to_eval, Matrix)==false
        nodes_to_eval = reshape( nodes_to_eval, (1,3) )
    end

    interp_vals = [ grid_int_interp( nodes_to_eval[i,1], nodes_to_eval[i,2], nodes_to_eval[i,3]) for i in (1:size(nodes_to_eval,1) ) ]

    interp_grad_vals = [ Interpolations.gradient( grid_int_interp, nodes_to_eval[i,1], nodes_to_eval[i,2], nodes_to_eval[i,3]) for i in (1:size(nodes_to_eval,1) ) ]

    return interp_vals, interp_grad_vals
end
