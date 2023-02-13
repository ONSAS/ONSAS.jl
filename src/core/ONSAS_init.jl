
using .StructuralAnalyses: StaticAnalysis, structure
using .StructuralSolvers: AbstractSolver

import .StructuralSolvers: init

"Returns the initialized analysis and solution struct. "
function init(sa::StaticAnalysis, alg::AbstractSolver, args...; kwargs...)

    s = structure(sa)

    λ0 = first(load_factors(sa))

    # Apply load BC into the global external forces vector 
    _apply_load_bc!(s, λ0)

    # Build initial external forces vector and apply boundary conditions
    _apply_disp_bc!(s, λ0)

    return sa
end