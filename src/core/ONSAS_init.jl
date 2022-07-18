"""
Function used to initialize the model solution and properties structs. 
"""
function ONSAS_init( materials, geometries, boundary_conditions, initial_conditions, mesh, analysis_settings )


    # ---------------------------------------
    # create properties
    neum_dofs = compute_neum_dofs( mesh, boundary_conditions )


    properties = ModelProperties( materials, geometries, boundary_conditions, neum_dofs, mesh, analysis_settings )
    # ---------------------------------------


    # ---------------------------------------
    # create initial solution
    num_nodes = size( mesh.nodal_coords, 1)

    U       = zeros( 6*num_nodes )
    Udot    = zeros( 6*num_nodes )
    Udotdot = zeros( 6*num_nodes ) # TO DO compute acceleration 
  
    # system_matrix, system_rhs = assemble_system( properties, Unp1k, neum_dofs  )
    # solution = ModelSolution( 0.0, U, Udot, Udotdot, system_matrix, system_rhs )
    solution = 1
    # ---------------------------------------


    return solution, properties
end





function generate_neum( nNodes )

    neumDofs = nodes2dofs( (1:nNodes), 6 )
    neumDofs[2:2:end] .= 0
    return neumDofs
end


function compute_neum_dofs( mesh, boundary_conditions)

    num_nodes = size( mesh.nodal_coords, 1)
    
    # keep non-zero boundary conds indexes
    boundary_types = sort( unique( mesh.MGBIValsMat[:,3]) )
    boundary_types[1] == 0 && popat!( boundary_types, 1)

    for BCnum in boundary_types

    end


    neum_dofs = generate_neum( num_nodes )


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

