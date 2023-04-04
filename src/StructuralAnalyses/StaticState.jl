using StaticArrays: @MVector
using SparseArrays: SparseMatrixCSC
using Dictionaries: Dictionary

using ...Meshes: num_dofs
using ...StructuralModel: AbstractStructure, num_free_dofs
using ...StructuralAnalyses: internal_forces, external_forces

import ..StructuralAnalyses: tangent_matrix, residual_forces
import ...StructuralSolvers: _solve, _update!, _reset!, external_forces

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