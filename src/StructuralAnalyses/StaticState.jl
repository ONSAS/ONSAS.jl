using SparseArrays: spzeros
using Dictionaries: Dictionary
using Reexport

using ...Meshes
using ...StructuralModel
using ...StructuralAnalyses
using ...StructuralSolvers
using ..Assemblers
using ...Utils
using ...Nodes

@reexport import ..StructuralAnalyses: tangent_matrix, residual_forces!, reset!

export StaticState

"""
Stores the relevant static variables of the structure during the displacements iteration.
"""
struct StaticState{DU<:AbstractVector,U<:AbstractVector,
                   FE<:AbstractVector,FI<:AbstractVector,K<:AbstractMatrix,
                   E<:Dictionary,S<:Dictionary} <: AbstractStructuralState
    "Free degrees of freedom."
    free_dofs::Vector{Dof}
    "Displacements vector increment."
    ΔUᵏ::DU
    "Displacements vector."
    Uᵏ::U
    "External forces vector."
    Fₑₓₜᵏ::FE
    "Internal forces vector."
    Fᵢₙₜᵏ::FI
    "System's tangent matrix."
    Kₛᵏ::K
    "Residual forces cache."
    res_forces::DU
    "Vector with straings for each element."
    ϵᵏ::E
    "Vector with stresses for each element."
    σᵏ::S
    "Assembler handler."
    assembler::Assembler
    "Current iteration state."
    iter_state::ResidualsIterationStep
    function StaticState(fdofs::Vector{Dof}, ΔUᵏ::DU, Uᵏ::U, Fₑₓₜᵏ::FE, Fᵢₙₜᵏ::FI, Kₛᵏ::K,
                         res_forces::DU, ϵᵏ::E,
                         σᵏ::S,
                         assembler::Assembler,
                         iter_state::ResidualsIterationStep) where {DU,U,FE,FI,K,E,S}
        # # Check dimensions
        @assert length(ΔUᵏ) == length(fdofs)
        @assert size(Kₛᵏ, 1) == size(Kₛᵏ, 2) == length(Fᵢₙₜᵏ) == length(Fₑₓₜᵏ) == length(Uᵏ)
        new{DU,U,FE,FI,K,E,S}(fdofs, ΔUᵏ, Uᵏ, Fₑₓₜᵏ, Fᵢₙₜᵏ, Kₛᵏ, res_forces, ϵᵏ, σᵏ, assembler,
                              iter_state)
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
    fdofs = free_dofs(s)
    return StaticState(fdofs, ΔUᵏ, Uᵏ, Fₑₓₜᵏ, Fᵢₙₜᵏ, Kₛᵏ, res_forces, ϵᵏ, σᵏ, assemblerᵏ,
                       iter_state)
end

function Base.show(io::IO, sc::StaticState)
    nu = length(sc.Uᵏ)
    K = sc.Kₛᵏ
    s = size(K)
    println("• StaticState with $nu-dofs displacements vector Uᵏ " *
            "and $(s[1]) × $(s[2]) tangent matrix Kₛᵏ with $(length(K.nzval)) stored entries.")
end

"Update and return the current residual forces of the `StaticState` `sc`."
function residual_forces!(sc::StaticState)
    return sc.res_forces .= view(external_forces(sc), free_dofs(sc)) -
                            view(internal_forces(sc), free_dofs(sc))
end

"Return the current system tangent matrix of the `StaticState` `sc`."
tangent_matrix(sc::StaticState) = sc.Kₛᵏ

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
