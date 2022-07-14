function ONSAS_init( materials, elements, boundary_conditions, initial_conditions, mesh, analysis_settings )


    solution = ModelSolution( 0.0 )

    properties = ModelProperties( mesh, analysis_settings )

    return solution, properties
end