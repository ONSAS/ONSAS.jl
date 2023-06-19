"""
Module defining solutions interface, `AbstractSolution`.

Each solution can at least access to displacements and external forces. More complete solutions
can include different forces vector, stresses and strains.

Of course the solution contains on the analysis and solver used to solve the problem.
"""
module Solutions

using Reexport

using ..Entities
using ..Meshes
using ..Interpolators
using ..Handlers
using ..Utils
using ..StructuralSolvers
using ..Nodes

@reexport import ..Entities: internal_forces, inertial_forces, strain, stress

export AbstractSolution, StatesSolution, stresses, strains, states, analysis, solver,
       displacements, external_forces, iteration_residuals

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
struct StatesSolution{ST<:Vector,A,SS<:AbstractSolver} <: AbstractSolution
    "Vector containing the converged structural states at each step."
    states::ST
    "Analysis solved."
    analysis::A
    "Solver employed."
    solver::SS
    "Constructor with empty `AbstractStructuralState`s `Vector` and type `S`."
    function StatesSolution(analysis::A, solver::SS) where {A,SS<:AbstractSolver}
        new{Vector{Any},A,SS}([], analysis, solver)
    end
end

"Show the states solution."
function Base.show(io::IO, ::MIME"text/plain", sol::StatesSolution)
    println("Analysis solved:")
    println("----------------\n")
    show(io, analysis(sol))

    println("\nSolver employed:")
    println("----------------\n")
    show(io, solver(sol))

    println("\nAccessors:")
    println("----------\n")
    println("• `displacements`")
    println("• `strain`")
    println("• `stress`")
    println("• `internal_forces`")
    println("• `external_forces`")
end
"Return the solved states."
states(sol::StatesSolution) = sol.states

for f in [:displacements, :internal_forces, :external_forces]
    "Return the $f vector Uᵏ at every time steps."
    @eval $f(st_sol::StatesSolution) = $f.(states(st_sol))

    "Return the $f at a certain dof for every time step."
    @eval $f(st_sol::StatesSolution, dof::Dof) = getindex.($f(st_sol), index(dof))

    "Return the $f at a certain dof's vector for every time step."
    @eval $f(st_sol::StatesSolution, vdof::Vector{Dof}) = [$f(st_sol, dof) for dof in vdof]

    "Return the $f at a certain node for every time step."
    @eval $f(st_sol::StatesSolution, n::AbstractNode) = $f(st_sol, reduce(vcat, collect(dofs(n))))

    "Return the $f component at a certain node for every time step."
    @eval $f(st_sol::StatesSolution, n::AbstractNode, component::Int) = $f(st_sol, n)[component]

    "Return the $f of an element for every time step."
    @eval $f(st_sol::StatesSolution, e::AbstractElement) = [$f(st_sol, n) for n in nodes(e)]
end

"Return the residuals iteration object at every time step."
iteration_residuals(st_sol::StatesSolution) = iteration_residuals.(states(st_sol))

for f in [:stress, :strain]
    "Return the $f for every time step."
    @eval $f(st_sol::StatesSolution) = $f.(states(st_sol))

    "Return the $f at a certain element for every time step."
    @eval function $f(st_sol::StatesSolution, e::AbstractElement)
        [getindex($f.(states(st_sol))[step], e) for step in 1:length(states(st_sol))]
    end
end

"Return the displacements solution at the points in the point evaluator handler."
function displacements(st_sol::StatesSolution, peh::PointEvalHandler)
    interpolate(st_sol, displacements, interpolator(peh))
end

"Return the displacements component in the point evaluator handler."
function displacements(st_sol::StatesSolution, peh::PointEvalHandler, i::Int)
    sol_points = getindex.(displacements(st_sol, peh), i)
    if length(sol_points) == 1
        getindex(sol_points)
    else
        sol_points
    end
end

"Return the internal forces solution at the points in the point eval handler."
function internal_forces(st_sol::StatesSolution, peh::PointEvalHandler)
    interpolate(st_sol, internal_forces, interpolator(peh))
end

"Return the internal force component at the points in the point eval handler."
function internal_forces(st_sol::StatesSolution, peh::PointEvalHandler, i::Int)
    getindex.(internal_forces(st_sol, peh), i)
end

"Return the external forces at the points in the point eval handler."
function external_forces(st_sol::StatesSolution, peh::PointEvalHandler)
    interpolate(st_sol, external_forces, interpolator(peh))
end

"Return the external force component at the points in the point eval handler."
function external_forces(st_sol::StatesSolution, peh::PointEvalHandler, i::Int)
    getindex.(external_forces(st_sol, peh), i)
end

"Return the stress at the elements where the points are located in the point eval handler."
function stress(st_sol::StatesSolution, peh::PointEvalHandler)
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
function strain(st_sol::StatesSolution, peh::PointEvalHandler)
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

end # module
