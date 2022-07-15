function ONSAS_init( materials, elements, boundary_conditions, initial_conditions, mesh, analysis_settings )

    welcome_message()

    num_nodes = size( mesh.nodal_coords, 1)

    U       = zeros( 6*num_nodes )
    Udot    = zeros( 6*num_nodes )
    Udotdot = zeros( 6*num_nodes ) # TO DO compute acceleration 
  
    solution = ModelSolution( 0.0, U, Udot, Udotdot )

    properties = ModelProperties( mesh, analysis_settings )

    return solution, properties
end



function welcome_message()
    print("\n =====================\n  WELCOME TO ONSAS.jl \n")
end