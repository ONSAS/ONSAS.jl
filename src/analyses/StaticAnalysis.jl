using StaticArrays: @MVector
using SparseArrays: SparseMatrixCSC

using ..Elements: local_dofs
using ..StructuralModel: AbstractStructure, num_dofs, num_elements
using ..StructuralSolvers: AbstractSolver, NewtonRaphson
using ..StructuralAnalyses: AbstractStructuralState, AbstractStructuralAnalysis
using ..StructuralAnalyses: Assembler

import ..Utils: _unwrap
import ..StructuralAnalyses: current_state, initial_time, current_time, final_time, next!
import ..StructuralSolvers: init, _solve

export StaticAnalysis, load_factors, current_load_factor

#==============#
# Static State
#==============#
"""
An `StaticState` object facilitates the process of storing the relevant static variables of the structure. 
### Fields:
- `ΔUᵏ`   -- stores displacements vector increment.
- `Uᵏ`    -- stores displacements vector.
- `Fₑₓₜᵏ` -- stores external forces vector.
- `Fᵢₙₜᵏ` -- stores internal forces vector.
- `Kₛᵏ`   -- stiffness tangent matrix of the structure. 
"""
struct StaticState <: AbstractStructuralState
    ΔUᵏ::AbstractVector
    Uᵏ::AbstractVector
    Fₑₓₜᵏ::AbstractVector
    Fᵢₙₜᵏ::AbstractVector
    Kₛᵏ::AbstractMatrix
    ϵ::AbstractVector
    σ::AbstractVector
    assembler::Assembler
end

"Returns a default static case for a given mesh."
function StaticState(s::AbstractStructure)
    n_dofs = num_dofs(s)
    n_elements = num_elements(s)
    Uᵏ = @MVector zeros(n_dofs)
    ΔUᵏ = similar(Uᵏ)
    Fₑₓₜᵏ = similar(Uᵏ)
    Fᵢₙₜᵏ = similar(Uᵏ)
    Kₛᵏ = SparseMatrixCSC(zeros(n_dofs, n_dofs))
    ϵᵏ = Vector{}(undef, n_elements)
    σᵏ = Vector{}(undef, n_elements)
    assemblerᵏ = Assembler(s)
    StaticState(ΔUᵏ, Uᵏ, Fₑₓₜᵏ, Fᵢₙₜᵏ, Kₛᵏ, ϵᵏ, σᵏ, assemblerᵏ)
end

residual_forces(sc::StaticState) = sc.Fₑₓₜᵏ - sc.Fᵢₙₜᵏ
systemΔu_matrix(sc::StaticState) = sc.Kₛᵏ

_unwrap(sc::StaticState) = (sc.ΔUᵏ, sc.Uᵏ, sc.Fₑₓₜᵏ, sc.Fᵢₙₜᵏ, sc.Kₛᵏ, sc.ϵᵏ, sc.σᵏ, sc.assemblerᵏ)

#================#
# Static Analysis
#================#
""" StaticAnalysis struct.
A `StaticAnalysis` is a collection of parameters for defining the static analysis of the structure. 

### Fields:
- `s`             -- Stores the structure to be analyzed.
- `state`         -- Stores the structural state.
- `λᵥ`            -- Stores the load factors vector of the analysis
- `current_step`  -- Stores the current load factor step
"""

struct StaticAnalysis <: AbstractStructuralAnalysis
    s::AbstractStructure
    state::StaticState
    λᵥ::Vector{<:Real}
    current_step::Int
end

function StaticAnalysis(s::AbstractStructure, λᵥ::Vector{<:Real}, current_step=0)
    StaticAnalysis(s, StaticState(s), λᵥ, current_step)
end

function StaticAnalysis(s::AbstractStructure, t₁::Real=1.0; NSTEPS=10, init_state=StaticState(s))
    t₀ = t₁ / NSTEPS
    λᵥ = LinRange(t₀, t₁, NSTEPS) |> collect
    StaticAnalysis(s, StaticState(s), λᵥ, 0)
end

initial_time(sa::StaticAnalysis) = first(load_factors(sa))

current_time(sa::StaticAnalysis) = sa.current_step

current_state(sa::StaticAnalysis) = sa.state

final_time(sa::StaticAnalysis) = last(load_factors(sa))

"Returns load factors vector"
load_factors(sa::StaticAnalysis) = sa.λᵥ

"Returns the current load factor"
current_load_factor(sa::StaticAnalysis) = load_factors(sa)[current_time(sa)]

function next!(sa::StaticAnalysis)
    next_step = current_load_factor(sa) + 1
    (next_step <= length(load_factors(sa))) || throw(ArgumentError("Analysis is done λᵏ = $(current_load_factor(sa))"))
    sa.current_step += 1
end

#================#
# Solve
#================#
"Returns the initialized analysis and solution struct. "
function init(sa::StaticAnalysis, alg::AbstractSolver, args...; kwargs...)

    s = structure(sa)

    λ0 = first(load_factors(sa))

    # Creates the external force vector
    _apply_load_bc!(s, sa, λ0)

    # Build initial external forces vector and apply boundary conditions
    _apply_disp_bc!(s, sa, λ0)

    return sa
end

"Internal function to solve different analysis problem"
function _solve(sa::StaticAnalysis, alg::AbstractSolver, args...; kwargs...)

    s = structure(sa)

    while !is_done(sa)

        # Computes system residual forces tangent system matrix    
        assemble_system!(s, sa, alg)

        while !is_step_converged(sa)

            # Increment U
            step!(sa, alg)

            # Re-computes system Δu and residual_forces
            assemble_system!(s, sa, alg)

            # Check convergence
            check_convergence!(sa, alg)

        end

        next!(sa)

    end


    return sol
end


"Computes system residual forces and tangent system matrix for the analysis"
function assemble_system!(s::AbstractStructure, sa::StaticAnalysis, ::NewtonRaphson)

    s = structure(sa)
    state = current_state(sa)

    for mat in structural_materials
        for e in elements(s)

            # Nodes and dofs of the element
            element_nodes_dofs = dofs(e)
            element_local_dofs = local_dofs(e)
            global_dofs = element_nodes_dofs[element_local_dofs]

            u_e = displacements(e, state)
            fᵢₙₜ_e, Kᵢₙₜ_e, σ_e, ϵ_e = internal_forces(mat, e, u_e)


            # Computes element residual forces tangent system matrix
            assemble!(e, sa)

            # Assembles element residual forces tangent system matrix into global residual forces tangent system matrix
            assemble!(s, e, sa)
        end
    end

end