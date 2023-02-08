using ..StructuralSolvers: AbstractSolver, ConvergenceSettings
using ..StructuralAnalyses: StaticAnalysis, structure
using ..StructuralModel: current_state, free_dofs_indexes, _unwrap

""" Newton-Raphson solver struct.

### Fields:
- `Δt`  -- Step size
- `tol` -- Numerical tolerances
"""

struct NewtonRaphson <: AbstractSolver
    Δt::Number
    tol::ConvergenceSettings
end

function step!(alg::NewtonRaphson, sa::StaticAnalysis, args...)

    # Current state
    s = structure(sa)
    Uᵏ, Fₑₓₜᵏ, Fᵢₙₜᵏ, Kₛᵏ = _unwrap(current_state(s))
    fdofs_indexes = free_dofs_indexes(s)

    # Solve system
    Kₛ_redᵏ = view(Kₛᵏ, fdofs_indexes, fdofs_indexes)
    Fₑₓₜ_redᵏ = view(Fₑₓₜᵏ, fdofs_indexes)
    Fᵢₙₜ_redᵏ = view(Fᵢₙₜᵏ, fdofs_indexes)

    ΔUᵏ = Kₛ_redᵏ \ (Fₑₓₜ_redᵏ - Fᵢₙₜ_redᵏ)

    # Computes Uk
    U_redᵏ = view(Uᵏ, fdofs_indexes) + ΔUᵏ
    Uᵏ[fdofs_indexes] = U_redᵏ

    return Uᵏ, ΔUᵏ
end