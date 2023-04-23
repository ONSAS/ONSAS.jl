using SparseArrays: spzeros
using Dictionaries: Dictionary

using ...Meshes: num_dofs
using ...StructuralModel: AbstractStructure, num_free_dofs
using ...StructuralAnalyses: displacements, Δ_displacements,
                             internal_forces, external_forces, iteration_residuals, stress, strain
using ...StructuralSolvers: _reset!

import ..StructuralAnalyses: tangent_matrix, residual_forces!, reset!
import ...StructuralSolvers: _update!

export StaticState

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
                   E<:Dictionary,S<:Dictionary} <: AbstractStructuralState
    # Structure
    s::ST
    #Displacements
    ΔUᵏ::DU
    Uᵏ::U
    #Forces
    Fₑₓₜᵏ::FE
    Fᵢₙₜᵏ::FI
    Kₛᵏ::K
    res_forces::DU
    ϵᵏ::E
    σᵏ::S
    # Iter
    assembler::Assembler
    iter_state::ResidualsIterationStep
    function StaticState(s::ST, ΔUᵏ::DU, Uᵏ::U, Fₑₓₜᵏ::FE, Fᵢₙₜᵏ::FI, Kₛᵏ::K, res_forces::DU, ϵᵏ::E,
                         σᵏ::S,
                         assembler::Assembler,
                         iter_state::ResidualsIterationStep) where {ST,DU,U,FE,FI,K,E,S}
        # # Check dimensions
        @assert length(ΔUᵏ) == num_free_dofs(s)
        @assert size(Kₛᵏ, 1) == size(Kₛᵏ, 2) == length(Fᵢₙₜᵏ) == length(Fₑₓₜᵏ) == length(Uᵏ) ==
                num_dofs(s)
        return new{ST,DU,U,FE,FI,K,E,S}(s, ΔUᵏ, Uᵏ, Fₑₓₜᵏ, Fᵢₙₜᵏ, Kₛᵏ, res_forces, ϵᵏ, σᵏ,
                                        assembler, iter_state)
    end
end

"Constructor for `StaticState` given an `AbstractStructure` `s`."
function StaticState(s::AbstractStructure,
                     iter_state::ResidualsIterationStep=ResidualsIterationStep())
    n_dofs = num_dofs(s)
    n_fdofs = num_free_dofs(s)
    Uᵏ = zeros(n_dofs)
    ΔUᵏ = zeros(n_fdofs)
    Fₑₓₜᵏ = zeros(n_dofs)
    Fᵢₙₜᵏ = zeros(n_dofs)
    Kₛᵏ = spzeros(n_dofs, n_dofs)
    res_forces = zeros(n_fdofs)
    # Initialize pairs strains 
    ϵᵏ = dictionary([Pair(e, Matrix{Float64}(undef, (3, 3))) for e in elements(s)])
    σᵏ = dictionary([Pair(e, Matrix{Float64}(undef, (3, 3))) for e in elements(s)])
    assemblerᵏ = Assembler(s)
    return StaticState(s, ΔUᵏ, Uᵏ, Fₑₓₜᵏ, Fᵢₙₜᵏ, Kₛᵏ, res_forces, ϵᵏ, σᵏ, assemblerᵏ, iter_state)
end

"Update and return the current residual forces of the `StaticState` `sc`."
function residual_forces!(sc::StaticState)
    return sc.res_forces .= view(external_forces(sc), free_dofs(sc)) -
                            view(internal_forces(sc), free_dofs(sc))
end

"Return the current system tangent matrix of the `StaticState` `sc`."
tangent_matrix(sc::StaticState) = sc.Kₛᵏ

"Updates displacements in the `StaticState` `sc` with a displacements increment vector `ΔU`."
function _update!(sc::StaticState, ΔU::AbstractVector)
    return sc.Uᵏ[free_dofs(sc)] .+= ΔU
end

"Resets the `StaticState` assembled magnitudes and the iteration state."
function reset!(state::StaticState)
    # Reset assembled magnitudes
    internal_forces(state) .= 0.0
    tangent_matrix(state)[findall(!iszero, tangent_matrix(state))] .= 0.0
    _reset!(assembler(state))
    # Reset the stress and strains dictionaries
    for (e, _) in pairs(stress(state))
        stress(state)[e] .= zeros(3, 3)
        strain(state)[e] .= zeros(3, 3)
    end
    # Reset ext force
    external_forces(state) .= 0.0
    # Reset iteration state
    displacements(state) .= 0.0
    Δ_displacements(state) .= 0.0
    _reset!(iteration_residuals(state))
    # Return state
    @info "The structural state has been reset."
    return state
end
