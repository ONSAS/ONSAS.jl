
function boundary_cond_processing(Materials, Geometries, Mesh, LoadsBC, DofsBC)

    Conec = Mesh.elem_nodal_connec # element connectivity matrix
    Nodes = Mesh.nodal_coords # nodes coordinate matrix
    nnodes = size(Nodes, 1)

    # Loads BC
    LoadsBCTypes = unique(Mesh.MGBI_mat[:, 3])
    LoadsBCTypes = deleteat!(LoadsBCTypes, findall(x -> x == 0, LoadsBCTypes))
    # Dofs BC
    DofsBCTypes = unique(Mesh.MGBI_mat[:, 4])
    DofsBCTypes = deleteat!(DofsBCTypes, findall(x -> x == 0, DofsBCTypes))
    # Element types
    elemTypes = unique(Mesh.MGBI_mat[:, 2])

    #
    factorLoadsFext = []
    loadFactorsFuncCell = []
    diri_dofs = []

    # loop over BC in Mesh
    for indBC = 1:length(LoadsBCTypes)
        println(indBC)
        BC = LoadsBCTypes[indBC]
        println(BC)
        mesh_with_BC = findall(x -> x == BC, Mesh.MGBI_mat[:, 3])
        println(mesh_with_BC)
        aux = elem2NodalLoads(Conec, Nodes, Mesh, Geometries, mesh_with_BC, LoadsBC[BC])
        push!(factorLoadsFext, aux)
        # push!(loadFactorsFuncCell, f(t)=LoadsBCTypes[BC].loadsTimeFactor * t)
    end
    println(factorLoadsFext)

    return factorLoadsFext
    # return Conec, Nodes, factorLoadsFextCell, loadFactorsFuncCell, diriDofs, neumDofs, KS, userLoadsFilename
end