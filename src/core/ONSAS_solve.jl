
using ..StructuralSolvers: NewtonRaphson
using ..StructuralAnalyses: StaticAnalysis

"Internal function to solve different analysis problem"
function _solve(sa::StaticAnalysis, alg::NewtonRaphson, args...; kwargs...)

    s = structure(sa)

    while !is_done(sa, λᵏ)

        # Increment external forces vector Fₑₓₜ  
        update_external_forces!(s, λᵏ)

        while !is_step_converged(sa)

            # Compute tangents and internal forces
            update_internal_forces!(s, sa)

            # Increment U
            step!(s, sa, alg)

            check_convergence!(sa, alg)

        end

        next!(sa)

    end


    return sol
end


