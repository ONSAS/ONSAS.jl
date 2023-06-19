module Solvers

using Reexport

using ..StructuralSolvers

@reexport import ..StructuralSolvers: tolerances

export NewtonRaphson

"""
Newton-Raphson solver struct.
"""
struct NewtonRaphson <: AbstractSolver
    "Numerical tolerances."
    tol::ConvergenceSettings
end

"Show a `NewtonRaphson` solver."
function Base.show(io::IO, nr::NewtonRaphson)
    println("Newton-Raphson solver with tolerances:")
    show(io, tolerances(nr))
end

"Constructor for `NewtonRaphson` solver with default tolerances."
NewtonRaphson() = NewtonRaphson(ConvergenceSettings())

"Return `NewtonRaphson` solver tolerances."
tolerances(nr::NewtonRaphson) = nr.tol

end # module
