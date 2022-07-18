"""
Function used to solve the problem and generate a vector of ModelSolution
"""
function ONSAS_solve( model_properties, initial_solution; verbosity=false )


    solutions = [ initial_solution ]


    curr_solution = initial_solution
    curr_time     = curr_solution.time

    final_time    = model_properties.analysis_settings.final_time
    #delta_time    = model_properties.analysis_settings.delta_time
    
    while final_time_reached( curr_solution.time, final_time ) == false

        verbosity && println( " curr_time: ", curr_time )

        next_solution = time_step_iteration( curr_solution, model_properties, verbosity=verbosity )

        push!( solutions, next_solution )

        # update
        curr_solution = next_solution
 
 #       verbosity && println( " OJOOO curr_time: ", curr_time, " final_time:", final_time )

    end

    return solutions
end



# check if final time was reached
final_time_reached( ct, ft ) = ct >= ft
