
using .StructuralAnalyses: StaticAnalysis, structure
using .StructuralSolvers: AbstractSolver

import .StructuralSolvers: init

"Returns load factors vector"
_load_factors_vector(sa::StaticAnalysis, alg::AbstractSolver) =
    LinRange(0, final_time(sa), ceil(Int, final_time(sa) / step_size(alg)))


"Returns the initialized analysis and solution struct. "
function init(sa::StaticAnalysis, alg::AbstractSolver, args...; kwargs...)

    s = structure(sa)
    #TODO : add an optinal load factors vector in kwargs 
    λs = _load_factors_vector(sa, alg)

    # Apply load BC into the global external forces vector 
    _apply_load_bc!(s, first(λs))

    # Build initial external forces vector and apply boundary conditions
    _apply_disp_bc!(s, first(λs))

    return sa
end