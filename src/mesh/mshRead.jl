


#joinpath("examples","cube.msh")

function mshRead( msh_filename )

    nodesCoordMat, connectivity, physicalNames, MEBIVec = MshFileReader( msh_filename )

    MEBIValsMat = parse_physicalNames( physicalNames )

    return nodesCoordMat, connectivity, MEBIValsMat, MEBIVec
end


function parse_physicalNames( physicalNames )

    print(physicalNames)

    return 1
end