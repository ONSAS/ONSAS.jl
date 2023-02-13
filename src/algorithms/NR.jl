using ..StructuralSolvers: AbstractSolver, ConvergenceSettings

export NewtonRaphson

""" Newton-Raphson solver struct.

### Fields:
- `tol` -- Numerical tolerances
"""

struct NewtonRaphson <: AbstractSolver
    tol::ConvergenceSettings
end

# Default tolerances
NewtonRaphson() = NewtonRaphson(ConvergenceSettings())