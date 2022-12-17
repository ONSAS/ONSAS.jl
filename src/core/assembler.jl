"""
assembler is the function in charge of the construction of the tangent matrix and the RHS vector of the main FEM iteration system 

The inputs are:
 - nodalCoords: a matrix that in the i-th colum contains the x, y and z coordinates as a column vector, for the i-th node of the mesh.
 - elemNodalConnec: a vector, that in the i-th entry contains a vector of the indexes of the node-connectivity of that element. 
 - MEBIValsMat: a matrix of indexes, at each column contains the the four MEBI parameters
 - 
"""
function assembler(materials, geometries, mesh, solution; verbosityBool=false)

    # General parameters
    elem_nodal_connec = mesh.elem_nodal_connec
    nodal_coords = mesh.nodal_coords

    num_elems = length(elem_nodal_connec)
    num_nodes = size(nodal_coords, 2)

    # Current solution
    isempty(currSol) ? currSol = zeros(6 * numNodes) : nothing

    # 
    Kmatrices = []
    # Sparse matrix indexes
    row_indexes = Vector{Int}()
    col_indexes = Vector{Int}()
    vals_vector = Vector{Float64}()

    # start assemblage
    for elem in 1:num_elems

        verbosityBool && println("element #: ", elem)

        # extract the MEBI entry
        elemMEBIEntry = view(MEBIVec, elem)

        elemMat = MEBIValsMat[elemMEBIEntry, 1]
        elemGeo = MEBIValsMat[elemMEBIEntry, 2]

        verbosityBool && println("  elemMebi entry", elemMEBIEntry)

        if elemMat > 0  # if it has a material

            nodesOfThisElem = elemNodalConnec[elem]
            dofsElem = nodes2dofs(nodesOfThisElem, 6)

            if cmp(geometries[elemEle].type, "truss") == 0

                dofsElem = dofsElem[[1, 5, 7, 11]]
                fintelem, Kelem = linear_truss(materialsData, elements_Geometry[elemEle], nodalCoords[nodesOfThisElem, :], currSol[dofsElem])

            elseif cmp(elementsData[elemEle].type, "tetrahedron") == 0

                dofsElem = dofsElem[Vector(1:2:6*4)]
                # E = materialsData[ MEBIValsMat[ elemMEBIEntry, 1 ] ].youngModulus
                Kelem = linear_tetrahedron(E)
            end

            elem_row_indexes, elem_col_indexes, elem_vals_vector = compute_sparse_indexes(Kelem, dofsElem)

            append!(row_indexes, elem_row_indexes)
            append!(col_indexes, elem_col_indexes)
            append!(vals_vector, elem_vals_vector)
            push!(Kmatrices, Kelem)
        end



    end

    KG = sparse(row_indexes, col_indexes, vals_vector, 6 * numNodes, 6 * numNodes)

    FG, neumDofs = assemble_fext(nodalCoords, elemNodalConnec, MEBIValsMat, MEBIVec, BCsData)

    return KG, FG, Kmatrices, neumDofs
end






