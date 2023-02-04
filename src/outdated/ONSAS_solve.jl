struct NewtonRaphson <: AbstractAlgorithm
    Δt::Float64
    final_time::Float64
end


"""
Function used to solve the problem and generate a vector of ModelSolution
"""
function ONSAS_solve(Model_properties, Initial_solution; verbosity=false)

    solutions = [Initial_solution]

    curr_solution = Initial_solution
    curr_time = curr_solution.time

    final_time = Model_properties.analysis_settings.final_time
    #delta_time    = Model_properties.analysis_settings.delta_time

    while final_time_reached(curr_solution.time, final_time) == false

        verbosity && println(" curr_time: ", curr_time)

        next_solution = time_step_iteration(curr_solution, Model_properties, verbosity=verbosity)

        push!(solutions, next_solution)

        # update
        curr_solution = next_solution

        #       verbosity && println( " OJOOO curr_time: ", curr_time, " final_time:", final_time )

    end

    return solutions
end



# check if final time was reached
final_time_reached(ct, ft) = ct >= ft
