module StaticAnalyses

using Dictionaries: dictionary
using Reexport: @reexport
using StaticArrays: @MVector
using SparseArrays: SparseMatrixCSC
using IterativeSolvers: cg
using LinearAlgebra: norm

@reexport using ..Materials
@reexport using ...Elements
@reexport using ...Meshes
@reexport using ...StructuralModel: AbstractStructure, num_dofs, num_elements, boundary_conditions, fixed_dof_bcs, load_bcs
@reexport using ..StructuralAnalyses
@reexport using ...StructuralSolvers

import ..StructuralAnalyses: _assemble!, initial_time, current_time, final_time, _next!,
    residual_forces, tangent_matrix, iteration_residuals, is_done, reset!

import ...StructuralSolvers: _solve, _update!, _step!, _reset!, external_forces

export StaticState
export StaticAnalysis, load_factors, current_load_factor

#==============#
# Static State
#==============#
"""
An `StaticState` object facilitates the process of storing the relevant static variables of the structure
during the displacements iteration. 
### Fields:
- `s`           -- stores the structure
- `ΔUᵏ`         -- stores displacements vector increment.
- `Uᵏ`          -- stores displacements vector.
- `Fₑₓₜᵏ`       -- stores external forces vector.
- `Fᵢₙₜᵏ`       -- stores internal forces vector.
- `Kₛᵏ`         -- stores the system tangent matrix.
- `ϵᵏ`          -- stores a vector with strains for each element.
- `σᵏ`          -- stores a vector with stresses for each element.
- `assembler`   -- assembler handler object 
- `iter_state`  -- current Δu iteration state 
"""
struct StaticState{ST<:AbstractStructure,
    DU<:AbstractVector,U<:AbstractVector,
    FE<:AbstractVector,FI<:AbstractVector,K<:AbstractMatrix,
    E<:Dictionary,S<:Dictionary
} <: AbstractStructuralState
    # Structure
    s::ST
    #Displacements
    ΔUᵏ::DU
    Uᵏ::U
    #Forces
    Fₑₓₜᵏ::FE
    Fᵢₙₜᵏ::FI
    Kₛᵏ::K
    ϵᵏ::E
    σᵏ::S
    # Iter
    assembler::Assembler
    iter_state::ResidualsIterationStep
    function StaticState(s::ST, ΔUᵏ::DU, Uᵏ::U, Fₑₓₜᵏ::FE, Fᵢₙₜᵏ::FI, Kₛᵏ::K, ϵᵏ::E, σᵏ::S,
        assembler::Assembler, iter_state::ResidualsIterationStep) where {ST,DU,U,FE,FI,K,E,S}
        # # Check dimensions
        @assert length(ΔUᵏ) == num_free_dofs(s)
        @assert size(Kₛᵏ, 1) == size(Kₛᵏ, 2) == length(Fᵢₙₜᵏ) == length(Fₑₓₜᵏ) == length(Uᵏ) == num_dofs(s)
        new{ST,DU,U,FE,FI,K,E,S}(s, ΔUᵏ, Uᵏ, Fₑₓₜᵏ, Fᵢₙₜᵏ, Kₛᵏ, ϵᵏ, σᵏ, assembler, iter_state)
    end
end

"Constructor for `StaticState` given an `AbstractStructure` `s`."
function StaticState(s::AbstractStructure)
    n_dofs = num_dofs(s)
    n_fdofs = num_free_dofs(s)
    Uᵏ = @MVector zeros(n_dofs)
    ΔUᵏ = @MVector zeros(n_fdofs)
    Fₑₓₜᵏ = @MVector zeros(n_dofs)
    Fᵢₙₜᵏ = similar(Fₑₓₜᵏ)
    Kₛᵏ = SparseMatrixCSC(zeros(n_dofs, n_dofs))
    # Initialize pairs strains 
    ϵᵏ = dictionary([Pair(e, Matrix{Float64}(undef, (3, 3))) for e in elements(s)])
    σᵏ = dictionary([Pair(e, Matrix{Float64}(undef, (3, 3))) for e in elements(s)])
    assemblerᵏ = Assembler(s)
    StaticState(s, ΔUᵏ, Uᵏ, Fₑₓₜᵏ, Fᵢₙₜᵏ, Kₛᵏ, ϵᵏ, σᵏ, assemblerᵏ, ResidualsIterationStep())
end

"Returns the current residual forces of the `StaticState` `sc`."
residual_forces(sc::StaticState) =
    external_forces(sc)[free_dofs(sc)] - internal_forces(sc)[free_dofs(sc)]

"Returns the current system tangent matrix of the `StaticState` `sc`."
tangent_matrix(sc::StaticState) = sc.Kₛᵏ

"Updates displacements in the `StaticState` `sc` with a displacements increment vector `ΔU`."
function _update!(sc::StaticState, ΔU::AbstractVector)
    sc.ΔUᵏ .= ΔU
    sc.Uᵏ[index.(free_dofs(sc))] .+= ΔU
end

"Resets the `StaticState` assembled magnitudes."
function _reset!(state::StaticState)
    _reset!(assembler(state))
    internal_forces(state) .= 0.0
    tangent_matrix(state)[findall(!iszero, tangent_matrix(state))] .= 0.0
end

#================#
# Static Analysis
#================#
""" StaticAnalysis struct.
A `StaticAnalysis` is a collection of parameters for defining the static analysis of the structure. 
In the static analysis, the structure is analyzed at a given load factor (this variable is analog to time).
### Fields:
- `s`             -- stores the structure to be analyzed.
- `state`         -- stores the structural state.
- `λᵥ`            -- stores the load factors vector of the analysis
- `current_step`  -- stores the current load factor step
"""
mutable struct StaticAnalysis{S<:AbstractStructure,LFV<:AbstractVector{<:Real}} <: AbstractStructuralAnalysis
    s::S
    state::StaticState
    λᵥ::LFV
    current_step::Int
    function StaticAnalysis(s::S, λᵥ::LFV; initial_step::Int=0) where {S<:AbstractStructure,LFV<:AbstractVector{<:Real}}
        new{S,LFV}(s, StaticState(s), λᵥ, initial_step)
    end
end

"Constructor for `StaticAnalysis` given a final time (or load factor) `t₁` and the number of steps `NSTEPS`."
function StaticAnalysis(s::AbstractStructure, t₁::Real=1.0; NSTEPS=10, initial_step::Int=1)
    t₀ = t₁ / NSTEPS
    λᵥ = LinRange(t₀, t₁, NSTEPS) |> collect
    StaticAnalysis(s, λᵥ, initial_step=initial_step)
end

"Returns the initial load factor of an `StaticAnalysis` `sa`."
initial_time(sa::StaticAnalysis) = first(load_factors(sa))

"Returns the current load factor of an `StaticAnalysis` `sa`."
current_time(sa::StaticAnalysis) = load_factors(sa)[sa.current_step]

"Returns the final load factor of an `StaticAnalysis` `sa`."
final_time(sa::StaticAnalysis) = last(load_factors(sa))

"Returns `true` if the `StaticAnalysis` `sa` is completed."
function is_done(sa::StaticAnalysis)
    is_done_bool = if sa.current_step > length(load_factors(sa))
        sa.current_step -= 1
        true
    else
        false
    end
end

"Returns the final load factor vector of an `StaticAnalysis` `sa`."
load_factors(sa::StaticAnalysis) = sa.λᵥ

"Returns the current load factor of an `StaticAnalysis` `sa`."
current_load_factor(sa::StaticAnalysis) = current_time(sa)

"Jumps to the next current load factor defined in the `StaticAnalysis` `sa`."
_next!(sa::StaticAnalysis) = sa.current_step += 1

"Sets the current load factor of the `StaticAnalysis` `sa` to the initial load factor."
reset!(sa::StaticAnalysis) = sa.current_step = 1

#================#
# Solve
#================#

"Solves an `StaticAnalysis` `sa` with an AbstractSolver `alg`."
function _solve(sa::StaticAnalysis, alg::AbstractSolver, args...; kwargs...)

    s = structure(sa)

    # Initialize solution
    sol = StatesSolution(sa, alg)

    # load factors iteration 
    while !is_done(sa)

        _reset!(current_iteration(sa))

        _apply!(sa, load_bcs(s)) # Compute Fext

        @debug external_forces(current_state(sa))

        while isconverged!(current_iteration(sa), tolerances(alg)) isa NotConvergedYet

            # Computes residual forces and tangent matrix    
            _assemble!(s, sa)

            @debug view(internal_forces(current_state(sa)), index.(free_dofs(s)))
            @debug residual_forces(current_state(sa))
            @debug tangent_matrix(current_state(sa))[index.(free_dofs(s)), index.(free_dofs(s))]

            # Increment U 
            _step!(sa, alg)

            @debug current_iteration(sa)
            @debug isconverged!(current_iteration(sa), tolerances(alg))

        end

        # Save current state
        push!(sol, current_state(sa))

        _next!(sa)

    end

    return sol
end

"Assembles the Structure `s` (internal forces) during the `StaticAnalysis` `sa`."
function _assemble!(s::AbstractStructure, sa::StaticAnalysis)

    state = current_state(sa)

    # Reset assembler
    _reset!(state)

    for (mat, mat_elements) in pairs(materials(s))
        for e in mat_elements

            # Global dofs of the element (dofs where K must be added)
            u_e = view(displacements(state), index.(local_dofs(e)))
            fᵢₙₜ_e, kₛ_e, σ_e, ϵ_e = internal_forces(mat, e, u_e)

            # Assembles the element internal magnitudes 
            _assemble!(state, fᵢₙₜ_e, e)
            _assemble!(state, kₛ_e, e)
            _assemble!(state, σ_e, ϵ_e, e)

        end

    end

    # Insert values in the assembler objet into the sysyem tangent stiffness matrix
    _end_assemble!(state)

end

"Computes ΔU for solving the `StaticAnalyses` `sa` with a `NewtonRaphson` method."
function _step!(sa::StaticAnalysis, ::NewtonRaphson)

    # Extract state info
    c_state = current_state(sa)
    f_dofs_indexes = index.(free_dofs(c_state))

    # Compute Δu
    r = residual_forces(c_state)
    K = view(tangent_matrix(c_state), f_dofs_indexes, f_dofs_indexes)
    ΔU = cg(K, r)

    # Compute norms
    norm_ΔU = norm(ΔU)
    rel_norm_ΔU = norm_ΔU / norm(displacements(c_state))
    norm_r = norm(r)
    rel_norm_r = norm_r / norm(external_forces(c_state))

    # Update displacements into the state
    _update!(c_state, ΔU)

    # Update iteration 
    _update!(current_iteration(sa), norm_ΔU, rel_norm_ΔU, norm_r, rel_norm_r)

end


"Pushes the current state `c_state` into the `StatesSolution` `st_sol`."
function Base.push!(st_sol::StatesSolution, c_state::StaticState)

    # Pointers
    s = structure(c_state)
    # Empty assembler since the info is stored in k
    assemblerᵏ = Assembler(s)

    # Deep copies 
    Uᵏ = deepcopy(displacements(c_state))
    ΔUᵏ = deepcopy(Δ_displacements(c_state))
    fₑₓₜᵏ = deepcopy(external_forces(c_state))
    fᵢₙₜᵏ = deepcopy(internal_forces(c_state))
    Kₛᵏ = deepcopy(tangent_matrix(c_state))
    σᵏ = dictionary([e => deepcopy(σ) for (e, σ) in pairs(stress(c_state))])
    ϵᵏ = dictionary([e => deepcopy(ϵ) for (e, ϵ) in pairs(strain(c_state))])
    iter_state = deepcopy(iteration_residuals(c_state))

    push!(states(st_sol), StaticState(s, ΔUᵏ, Uᵏ, fₑₓₜᵏ, fᵢₙₜᵏ, Kₛᵏ, ϵᵏ, σᵏ, assemblerᵏ, iter_state))
end



end # module