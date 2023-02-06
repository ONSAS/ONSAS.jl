
# ================================
# Numerical convergence parameters
# ================================




struct NewtonRaphson <: AbstractAlgorithm
    tolk::Float64
    tolu::Float64
    tolf::Float64
    loadFactors::Vector{Float64}
    nTimes::Int64
    function NewtonRaphson(tolk, tolu, tolf, loadFactors)
        nTimes = length(loadFactors)
        return new(tolk, tolu, tolf, loadFactors, nTimes)
    end
end




# ================================
# Numerical convergence parameters
# ================================

const DEFAULT_RELATIVE_DISP_TOLERANCE = 1e-6
const DEFAULT_RELATIVE_FORCE_TOLERANCE = 1e-6
const DEFAULT_MAX_ITER = 20

""" ConvergenceSettings struct.

Facilitates the process of defining and checking numerical convergence. 
### Fields:

- `rel_disp_tol`    -- Relative displacement tolerance.
- `rel_force_tol`   -- Relative residual force tolerance.
- `max_iter`        -- Maximum number of iterations.
"""
Base.@kwdef struct ConvergenceSettings

    rel_disp_tol::Number = DEFAULT_RELATIVE_DISP_TOLERANCE
    rel_force_tol::Number = DEFAULT_RELATIVE_FORCE_TOLERANCE
    max_iter::Integer = DEFAULT_MAX_ITER

end