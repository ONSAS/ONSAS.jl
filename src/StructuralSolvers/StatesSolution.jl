using ..StructuralSolvers: AbstractSolution
using ..Meshes: PointEvalHandler, interpolator, points

import ..Elements: internal_forces, inertial_forces, strain, stress

export StatesSolution, stresses, strains, states, analysis, solver, iteration_residuals

struct StatesSolution{ST<:Vector,A,SS<:AbstractSolver} <: AbstractSolution
    states::ST
    analysis::A
    solver::SS
    "Constructor with empty `AbstractStructuralState`s `Vector` and type `S`."
    function StatesSolution(analysis::A, solver::SS) where {A,SS<:AbstractSolver}
        new{Vector{Any},A,SS}([], analysis, solver)
    end
end


"Returns the solved `AbstractStrcturalState`s. "
states(sol::StatesSolution) = sol.states

"Returns the `AbstractAnalysis` solved. "
analysis(sol::StatesSolution) = sol.analysis

"Returns the `AbstractSolver` solved. "
solver(sol::StatesSolution) = sol.solver


for f in [:displacements, :internal_forces, :external_forces]
    "Returns the $f vector Uᵏ at every time step."
    @eval $f(st_sol::StatesSolution) = $f.(states(st_sol))

    "Returns the $f of the `Dof` at every time step."
    @eval $f(st_sol::StatesSolution, dof::Dof) = getindex.($f(st_sol), index(dof))

    "Returns the a $f `Vector` at a `Vector` of `Dof`s at every time step."
    @eval $f(st_sol::StatesSolution, vdof::Vector{Dof}) = [$f(st_sol, dof) for dof in vdof]

    "Returns the $f of a `Node` `n` every time step."
    @eval $f(st_sol::StatesSolution, n::AbstractNode) = $f(st_sol, reduce(vcat, collect(dofs(n))))

    "Returns the $f of a `Element` `e` at every time step."
    @eval $f(st_sol::StatesSolution, e::AbstractElement) = [$f(st_sol, n) for n in nodes(e)]
end

"Returns the `IterationResidual`s object at every time step."
iteration_residuals(st_sol::StatesSolution) = iteration_residuals.(states(st_sol))

for f in [:stress, :strain]
    "Returns the $f at every time step."
    @eval $f(st_sol::StatesSolution) = $f.(states(st_sol))

    "Returns the $f of an `Element` `e` every time step."
    @eval $f(st_sol::StatesSolution, e::AbstractElement) = [getindex($f.(states(st_sol))[step], e) for step in 1:length(states(st_sol))]

end

"Returns the displacements solution  at the `PointEvalHandler` `peh`."
function displacements(st_sol::StatesSolution, peh::PointEvalHandler)

    interpolators = interpolator(peh)
    vec_points = points(peh)
    num_points = length(vec_points)

    sol_points = Vector{Vector{Vector{Float64}}}(undef, num_points)

    for index_p in 1:num_points
        interpolator_p = interpolators[index_p]

        # Compute the first node contribution and the sum up
        node_values = [weight * displacements(st_sol, node) for (node, weight) in pairs(interpolator_p)]
        p_values = reduce(+, node_values)
        sol_points[index_p] = p_values
    end
    return sol_points
end

"Returns the internal forces solution  at the `PointEvalHandler` `peh`."
function internal_forces(st_sol::StatesSolution, peh::PointEvalHandler)

    interpolators = interpolator(peh)
    vec_points = points(peh)
    num_points = length(vec_points)

    sol_points = Vector{Vector{Vector{Float64}}}(undef, num_points)

    for index_p in 1:num_points
        interpolator_p = interpolators[index_p]

        # Compute the first node contribution and the sum up
        node_values = [weight * internal_forces(st_sol, node) for (node, weight) in pairs(interpolator_p)]
        p_values = reduce(+, node_values)
        sol_points[index_p] = p_values
    end
    return sol_points
end

"Returns the external forces solution  at the `PointEvalHandler` `peh`."
function external_forces(st_sol::StatesSolution, peh::PointEvalHandler)

    interpolators = interpolator(peh)
    vec_points = points(peh)
    num_points = length(vec_points)

    sol_points = Vector{Vector{Vector{Float64}}}(undef, num_points)

    for index_p in 1:num_points
        interpolator_p = interpolators[index_p]

        # Compute the first node contribution and the sum up
        node_values = [weight * external_forces(st_sol, node) for (node, weight) in pairs(interpolator_p)]
        p_values = reduce(+, node_values)
        sol_points[index_p] = p_values
    end
    return sol_points
end