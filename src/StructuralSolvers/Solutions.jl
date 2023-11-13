"""
Module defining solutions interface, `AbstractSolution`.

Each solution can at least access to displacements and external forces. More complete solutions
can include different forces vector, stresses and strains.

Of course the solution contains on the analysis and solver used to solve the problem.
"""
module Solutions

using Reexport, PrettyTables, Dictionaries

using ..Utils
using ..Entities
using ..Nodes
using ..Meshes
using ..Interpolators
using ..Handlers
using ..Meshes
using ..Structures
using ..StructuralAnalyses
using ..StaticStates
using ..StructuralSolvers

@reexport import ..StructuralAnalyses: displacements, external_forces, iteration_residuals
@reexport import ..StructuralSolvers: residual_forces_tol, displacement_tol, criterion, iterations
@reexport import ..Entities: internal_forces, inertial_forces, strain, stress

export AbstractSolution, Solution, stresses, strains, states, analysis, solver,
       displacements, external_forces, iteration_residuals, deformed_node_positions

"""
Abstract supertype for all structural analysis solutions.

**Abstract Methods**
* [`displacements`](@ref)
* [`external_forces`](@ref)
* [`internal_forces`](@ref)
* [`stress`](@ref)
* [`strain`](@ref)

**Abstract fields**
* analysis
* solver
"""
abstract type AbstractSolution end

"Return the analysis solved."
analysis(sol::AbstractSolution) = sol.analysis

"Return the solver used to solve the analysis."
solver(sol::AbstractSolution) = sol.solver

"""
Solution that stores all intermediate arrays during the analysis.
"""
struct Solution{ST<:AbstractStaticState,A<:AbstractStructuralAnalysis,
                SS<:Union{AbstractSolver,Nothing}} <: AbstractSolution
    "Vector containing the converged structural states at each step."
    states::Vector{ST}
    "Analysis solved."
    analysis::A
    "Solver employed."
    solver::SS
end

"Constructor with empty `AbstractStructuralState`s `Vector` and type `S`."
function Solution(analysis::A,
                  solver::SS) where {A<:AbstractStructuralAnalysis,
                                     SS<:Union{AbstractSolver,Nothing}}

    # TODO Use concrete types.
    state = current_state(analysis)
    Uᵏ = displacements(state)
    ϵᵏ = strain(state)
    σᵏ = stress(state)

    #TODO: Create a general way to pre-allocate with the number of analysis steps
    states = Vector{StaticState}(undef, length(analysis.λᵥ))

    for i in 1:length(analysis.λᵥ)
        sol_Uᵏ = similar(Uᵏ)
        sol_σᵏ = dictionary([e => similar(σ) for (e, σ) in pairs(σᵏ)])
        sol_ϵᵏ = dictionary([e => similar(ϵ) for (e, ϵ) in pairs(ϵᵏ)])
        states[i] = StaticState(sol_Uᵏ, sol_ϵᵏ, sol_σᵏ)
    end

    Solution{StaticState,A,SS}(states, analysis, solver)
end

"Show the states solution."
function Base.show(io::IO, ::MIME"text/plain", solution::Solution)
    println("Analysis solved:")
    println("----------------\n")
    show(io, analysis(solution))

    println("\nSolver employed:")
    println("----------------\n")
    show(io, solver(solution))

    println("\nStats:")
    println("----------")
    # Check convergence
    is_any_step_not_converged = any([criterion_step isa Union{NotConvergedYet,MaxIterCriterion}
                                     for criterion_step in criterion(solution)])

    num_iterations = reduce(+, iterations(solution))
    avg_iterations = round(num_iterations / length(states(solution)); digits=1)
    println("• Number of linear systems solved: $num_iterations")
    println("• Average of iterations per step : $avg_iterations")
    println("• Convergence success            : $(!is_any_step_not_converged)")
    _print_table(solution)
end

"Print the solution table"
function _print_table(solution::AbstractSolution)
    header = ["iter", "time", "||Uᵏ||", "||ΔUᵏ||/||Uᵏ||", "||ΔRᵏ||", "||ΔRᵏ||/||Fₑₓₜ||",
              "convergence criterion", "iterations"]

    ΔU_rel = getindex.(displacement_tol(solution), 1)
    ΔU = getindex.(displacement_tol(solution), 2)
    ΔR_rel = getindex.(residual_forces_tol(solution), 1)
    ΔR = getindex.(residual_forces_tol(solution), 2)
    t = analysis(solution).λᵥ
    criterions = [string(criterion_step)[1:(end - 2)] for criterion_step in criterion(solution)]
    iters = iterations(solution)
    num_times = length(t)
    data = hcat(collect(1:num_times), t, ΔU, ΔU_rel, ΔR, ΔR_rel, criterions, iters)

    hl = Highlighter(;
                     f=(data, i, j) -> i % 2 == 0,
                     crayon=Crayon(; foreground=:white, background=:black, bold=:true))

    pretty_table(data;
                 highlighters=hl,
                 header=header,
                 display_size=(20, 1000),
                 vcrop_mode=:middle,
                 formatters=(ft_printf("%i", [1]),
                             ft_printf("%.2f", [2]),
                             ft_printf("%e", [3, 4, 5, 6])),
                 alignment=:c)
end

"Return the solved states."
states(sol::Solution) = sol.states

for f in [:displacements, :internal_forces, :external_forces]
    "Return the $f vector Uᵏ at every time steps."
    @eval $f(st_sol::Solution) = $f.(states(st_sol))

    "Return the $f at a certain dof for every time step."
    @eval $f(st_sol::Solution, dof::Dof) = getindex.($f(st_sol), Utils.index(dof))

    "Return the $f at a certain dof's vector for every time step."
    @eval $f(st_sol::Solution, vdof::Vector{Dof}) = [$f(st_sol, dof) for dof in vdof]

    "Return the $f at a certain node for every time step."
    @eval $f(st_sol::Solution, n::AbstractNode) = $f(st_sol, reduce(vcat, collect(dofs(n))))

    "Return the $f component at a certain node for every time step."
    @eval $f(st_sol::Solution, n::AbstractNode, component::Int) = $f(st_sol, n)[component]

    "Return the $f of an element for every time step."
    @eval $f(st_sol::Solution, e::AbstractElement) = [$f(st_sol, n) for n in nodes(e)]
end

for f in [:stress, :strain]
    "Return the $f for every time step."
    @eval $f(st_sol::Solution) = $f.(states(st_sol))

    "Return the $f at a certain element for every time step."
    @eval function $f(st_sol::Solution, e::AbstractElement)
        [getindex($f.(states(st_sol))[step], e) for step in 1:length(states(st_sol))]
    end
end

"Return the residuals iteration object at every time step."
iteration_residuals(st_sol::Solution) = iteration_residuals.(states(st_sol))

"Return the final residuals forces at each time step."
residual_forces_tol(st_sol::Solution) = residual_forces_tol.(iteration_residuals(st_sol))

"Return the final residual displacement at each time step."
displacement_tol(st_sol::Solution) = displacement_tol.(iteration_residuals(st_sol))

"Return the convergence criterion at each time step."
criterion(st_sol::Solution) = criterion.(iteration_residuals(st_sol))

"Return the number of iterations at each time step."
iterations(st_sol::Solution) = iterations.(iteration_residuals(st_sol))

"Return the displacements solution at the points in the point evaluator handler."
function displacements(st_sol::Solution, peh::PointEvalHandler)
    interpolate(st_sol, displacements, interpolator(peh))
end

"Return the displacements component in the point evaluator handler."
function displacements(st_sol::Solution, peh::PointEvalHandler, i::Int)
    sol_points = getindex.(displacements(st_sol, peh), i)
    if length(sol_points) == 1
        getindex(sol_points)
    else
        sol_points
    end
end

"Return the internal forces solution at the points in the point eval handler."
function internal_forces(st_sol::Solution, peh::PointEvalHandler)
    interpolate(st_sol, internal_forces, interpolator(peh))
end

"Return the internal force component at the points in the point eval handler."
function internal_forces(st_sol::Solution, peh::PointEvalHandler, i::Int)
    getindex.(internal_forces(st_sol, peh), i)
end

"Return the external forces at the points in the point eval handler."
function external_forces(st_sol::Solution, peh::PointEvalHandler)
    interpolate(st_sol, external_forces, interpolator(peh))
end

"Return the external force component at the points in the point eval handler."
function external_forces(st_sol::Solution, peh::PointEvalHandler, i::Int)
    getindex.(external_forces(st_sol, peh), i)
end

"Return the stress at the elements where the points are located in the point eval handler."
function stress(st_sol::Solution, peh::PointEvalHandler)
    point_to_element_vec = points_to_element(interpolator(peh))
    vec_points = points(peh)
    num_points = length(vec_points)

    sol_points = Vector{Vector{<:AbstractMatrix{Float64}}}(undef, num_points)

    for index_p in 1:num_points
        element_p = point_to_element_vec[index_p]
        sol_points[index_p] = stress(st_sol, element_p)
    end
    sol_points
end

"Return the strains at the elements where the points are located in the point eval handler."
function strain(st_sol::Solution, peh::PointEvalHandler)
    point_to_element_vec = points_to_element(interpolator(peh))
    vec_points = points(peh)
    num_points = length(vec_points)

    sol_points = Vector{Vector{<:AbstractMatrix{Float64}}}(undef, num_points)

    for index_p in 1:num_points
        element_p = point_to_element_vec[index_p]
        sol_points[index_p] = strain(st_sol, element_p)
    end
    sol_points
end

"Return deformed nodes of the mesh at time index."
function deformed_node_positions(st_sol::Solution, t_i::Int)
    mesh_nodes = nodes(mesh(analysis(st_sol).s))
    # Get nodes coordinates type
    CT = typeof(coordinates(rand(mesh_nodes)))

    deformed_points = Vector{CT}(undef, length(mesh_nodes))
    u_node = Vector{Float64}(undef, dimension(first(mesh_nodes)))

    for (i, node) in enumerate(mesh_nodes)
        displacements_node = displacements(st_sol, node)
        for dim in 1:dimension(node)
            u_node[dim] = displacements_node[dim][t_i]
        end
        deformed_points[i] = coordinates(node) + u_node
    end
    deformed_points
end

end # module
