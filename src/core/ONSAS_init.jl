"""
Function used to initialize the model solution and properties structs. 
"""
function ONSAS_init(Materials, Geometries, LoadsBC, DofsBC, Initial_conditions, Mesh, ConvergenceSettings, Algorithm)

    # ---------------------------------------
    # create properties
    # neum_dofs = compute_neum_dofs(Mesh, Boundary_conditions)
    neum_dofs = []
    Properties = ModelProperties(Materials, Geometries, LoadsBC, DofsBC, Mesh, ConvergenceSettings, Algorithm, neum_dofs)
    # ---------------------------------------

    # ---------------------------------------
    # boundary conds
    Conec = boundary_cond_processing(Materials, Geometries, Mesh, LoadsBC, DofsBC)
    # ---------------------------------------

    # ---------------------------------------
    # create initial solution
    nnodes = size(Mesh.nodal_coords, 1)

    U = zeros(6 * nnodes)
    U̇ = zeros(6 * nnodes)
    Ü = zeros(6 * nnodes) # TO DO compute acceleration 

    # system_matrix, system_rhs = assemble_system(Properties, Unp1k, neum_dofs)
    system_matrix, system_rhs = assemble_system(Properties, U, neum_dofs)
    # Solution = ModelSolution( 0.0, U, Udot, Udotdot, system_matrix, system_rhs )
    Solution = 1
    # ---------------------------------------

    return Solution, Properties
end





function generate_neum(nNodes)

    neumDofs = nodes2dofs(collect(1:nNodes), 6)
    neumDofs[2:2:end] .= 0
    return neumDofs
end


function compute_neum_dofs(Mesh, Boundary_conditions)

    nnodes = size(Mesh.nodal_coords, 1)

    # keep non-zero boundary conds indexes
    boundary_types = sort(unique(Mesh.MGBI_mat[:, 3]))
    boundary_types[1] == 0 && popat!(boundary_types, 1)

    for BCnum in boundary_types

    end


    neum_dofs = generate_neum(nnodes)


    # print("neumdofs", neum_dofs)
    # FG = zeros( 6*numNodes )

    # elementsBC = computeElemsByMEBI( MEBIValsMat[:,3], MEBIVec )

    # for indBC in (1:length(BCsData))
    #     print("indBC", indBC)
    #     for elem in elementsBC[indBC]
    #         print("elem:", elem)
    #         nodesElem = elemNodalConnec[ elem ]
    #         dofs      = nodes2dofs( nodesElem, 6)

    #         if length( BCsData[indBC].NeumannNodalDOFs) >0
    #             dofsNeu = dofs[ BCsData[indBC].NeumannNodalDOFs ]
    #             FG[ dofsNeu ] = FG[ dofsNeu ] + BCsData[indBC].NeumannNodalVals
    #         end

    #         if length( BCsData[indBC].DirichletNodalDOFs) >0
    #             @assert sum( BCsData[indBC].DirichletNodalVals )==0 "only homogeneous dirichlet BC are allowed by now!"
    #             print("dofs:", dofs)
    #             neumDofs[ dofs[ BCsData[indBC].DirichletNodalDOFs ] ] .= 0
    #         end
    #     end
    # end

    # print("neumdofs", neumDofs)

    # neumDofs = sort( unique(neumDofs) )
    # neumDofs[1] == 0 && popfirst!(neumDofs)

    # print("neumdofs", neumDofs)
    # print("listo\n")

    return neum_dofs
end

