
function boundary_cond_processing(Mesh, Materials, elements, boundaryConds, analysisSettings)

    Conec = Mesh.elem_nodal_connec # element connectivity matrix
    Nodes = Mesh.nodal_coords # nodes coordinate matrix
    nnodes = size(Nodes, 1)

    LoadsBCTypes = unique(Mesh.MGBI_mat[:, 3])
    SupportsBCTypes = unique(Mesh.MGBI_mat[:, 4])
    elemTypes = unique(Mesh.MGBI_mat[:, 2])

    LoadsBCTypes = deleteat!(LoadsBCTypes, findall(x -> x == 0, LoadsBCTypes))
    SupportsBCTypes = deleteat!(SupportsBCTypes, findall(x -> x == 0, SupportsBCTypes))

    factorLoadsFext = []
    loadFactorsFuncCell = []
    diri_dofs = []

    # loop over BC in Mesh
    for indBC = 1:length(LoadsBCTypes)
        BC = LoadsBCTypes[indBC]
        aux = elem2NodalLoads(Mesh, Nodes, boundaryConds[BC], BC, elements)
        push!(factorLoadsFext, aux)
        push!(loadFactorsFuncCell, f(t)=LoadsBCTypes[BC].loadsTimeFactor * t)
    end

    return Conec, Nodes, factorLoadsFextCell, loadFactorsFuncCell, diriDofs, neumDofs, KS, userLoadsFilename
end