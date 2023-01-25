function elem2NodalLoads(Conec, Nodes, Mesh, Geometries, mesh_with_BC, BC)

    nnodes = size(Nodes, 1)
    Fext = zeros(6 * nnodes)

    loadVals = BC.loadsBaseVals
    loadsCoordSystem = BC.loadsCoordSystem

    loadedNodes = Vector{Vector{Float64}}()

    for i in 1:length(mesh_with_BC)

        println(i)
        row = mesh_with_BC[i]
        println(row)
        elem_index = Mesh.MGBI_mat[row, 2]
        println(elem_index)
        elem_type = Geometries[elem_index].type
        println(elem_type)
        loadsVector = generate_loads(elem_type::AbstractElement, loadVals, loadsCoordSystem, Conec, Nodes, elem_index, Geometries)
        push!(loadedNodes, loadsVector)
    end

    for i in 1:size(loadedNodes, 1)
        dofs = nodes2dofs(loadedNodes[i][1], 6)
        Fext[dofs] = Fext[dofs] + loadedNodes[i][2:7]
    end
    return Fext
end

function generate_loads(elem_type::Node, loadVals, loadCoordSystem, Conec, Nodes, elem_index, Geometries)

    if cmp(loadCoordSystem, "Global") == 0
        nodes = Conec[elem_index]
    else
        error("Not implemented yet")
    end
    println(nodes)
    return vcat(nodes, loadVals)
end