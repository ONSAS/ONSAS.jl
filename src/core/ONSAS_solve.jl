
using ..StructuralSolvers: NewtonRaphson
using ..StructuralAnalyses: StaticAnalysis

"Internal function to solve different analysis problem"
function _solve(sa::StaticAnalysis, alg::NewtonRaphson, args...; kwargs...)

    # initial conditions
    u0 = prob.u0
    tspan = prob.tspan

    # solve
    sol = solve(prob, method, options)

    # postprocess
    postprocess(prob, sol, options)

    return sol
end


