"""
Module defining static states.
Each static state defines the structural state of an structure into a StaticAnalysis. This
state contains external and internal forces, displacements, stresses and strains.
"""
module StaticStates

using SparseArrays
using Dictionaries
using InteractiveUtils
using LinearAlgebra
using Reexport
using LinearSolve

using ..Entities
using ..Meshes
using ..Structures
using ..StructuralAnalyses
using ..StructuralSolvers
using ..Assemblers
using ..Utils
using ..Nodes

@reexport import ..Assemblers: reset!
@reexport import ..StructuralAnalyses: residual_forces!
@reexport import ..StructuralSolvers: tangent_matrix

export FullStaticState, StaticState

"""
Stores the relevant static variables of the structure during the displacements iteration.
"""
struct FullStaticState{DU<:AbstractVector,U<:AbstractVector,
                       FE<:AbstractVector,FI<:AbstractVector,K<:AbstractMatrix,
                       E<:Dictionary,S<:Dictionary} <: AbstractStaticState
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
    "Vector with strains for each element."
    ϵᵏ::E
    "Vector with stresses for each element."
    σᵏ::S
    "Assembler handler."
    assembler::Assembler
    "Current iteration state."
    iter_state::ResidualsIterationStep
    "Linear system cache"
    linear_system::LinearSolve.LinearCache
    function FullStaticState(fdofs::Vector{Dof},
                             ΔUᵏ::DU, Uᵏ::U,
                             Fₑₓₜᵏ::FE, Fᵢₙₜᵏ::FI,
                             Kₛᵏ::K, res_forces::DU,
                             ϵᵏ::E, σᵏ::S,
                             assembler::Assembler,
                             iter_state::ResidualsIterationStep,
                             linear_system::LinearSolve.LinearCache) where {DU,U,FE,FI,K,E,S}
        # Check dimensions
        @assert length(ΔUᵏ) == length(fdofs) == length(res_forces)
        @assert size(Kₛᵏ, 1) == size(Kₛᵏ, 2) == length(Fᵢₙₜᵏ) == length(Fₑₓₜᵏ) == length(Uᵏ)
        # Initialize linear system K.ΔU = R
        new{DU,U,FE,FI,K,E,S}(fdofs, ΔUᵏ, Uᵏ, Fₑₓₜᵏ, Fᵢₙₜᵏ, Kₛᵏ, res_forces, ϵᵏ, σᵏ, assembler,
                              iter_state, linear_system)
    end
end

"Default constructor for static state given an structure and iteration state."
function FullStaticState(s::AbstractStructure,
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
    ϵᵏ = dictionary([Pair(e, Symmetric(Matrix{Float64}(undef, (3, 3)))) for e in elements(s)])
    σᵏ = dictionary([Pair(e, Matrix{Float64}(undef, (3, 3))) for e in elements(s)])
    cache = dictionary(nameof(T) => elements_cache(T) for T in subtypes(AbstractElement))
    assemblerᵏ = Assembler(s, cache)
    fdofs = free_dofs(s)
    linear_system = init(LinearProblem(Kₛᵏ[fdofs, fdofs], res_forces))
    FullStaticState(fdofs, ΔUᵏ, Uᵏ, Fₑₓₜᵏ, Fᵢₙₜᵏ, Kₛᵏ, res_forces, ϵᵏ, σᵏ, assemblerᵏ, iter_state,
                    linear_system)
end

function Base.show(io::IO, sc::FullStaticState)
    nu = length(sc.Uᵏ)
    K = sc.Kₛᵏ
    s = size(K)
    println("• FullStaticState with $nu-dofs displacements vector Uᵏ " *
            "and $(s[1]) × $(s[2]) tangent matrix Kₛᵏ with $(length(K.nzval)) stored entries.")
end

"Update and return the current residual forces of the static state."
function residual_forces!(sc::FullStaticState)
    return sc.res_forces .= view(external_forces(sc), free_dofs(sc)) -
                            view(internal_forces(sc), free_dofs(sc))
end

"Return the current system tangent matrix form the static state ."
tangent_matrix(sc::FullStaticState) = sc.Kₛᵏ

"Reset the static state assembled magnitudes and the iteration state."
function reset!(state::FullStaticState)
    # Reset assembled magnitudes
    internal_forces(state) .= 0.0
    tangent_matrix(state)[findall(!iszero, tangent_matrix(state))] .= 0.0
    reset!(assembler(state))
    # Reset the stress and strains dictionaries
    for (e, _) in pairs(stress(state))
        stress(state)[e] .= Symmetric(zeros(3, 3))
        strain(state)[e] .= Symmetric(zeros(3, 3))
    end
    # Reset ext force
    external_forces(state) .= 0.0
    # Reset iteration state
    displacements(state) .= 0.0
    Δ_displacements(state) .= 0.0
    reset!(iteration_residuals(state))
    # Return state
    @info "The structural state has been reset."
    state
end

struct StaticState{U<:AbstractVector,E<:Dictionary,S<:Dictionary} <: AbstractStaticState
    "Displacements vector."
    Uᵏ::U
    "Vector with strains for each element."
    ϵᵏ::E
    "Vector with stresses for each element."
    σᵏ::S
    function StaticState(Uᵏ::U, ϵᵏ::E, σᵏ::S) where {U,E,S}
        new{U,E,S}(Uᵏ, ϵᵏ, σᵏ)
    end
end

function Base.show(io::IO, sc::StaticState)
    nu = length(sc.Uᵏ)
    println("• StaticState with $nu-dofs displacements vector Uᵏ.")
end

end # module
