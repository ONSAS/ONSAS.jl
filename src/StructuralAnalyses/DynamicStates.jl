"""
Module defining Dynamic states.
Each Dynamic state defines the structural state of an structure into a DynamicAnalysis. This
state contains external and internal forces, displacements, stresses and strains.
"""
module DynamicStates

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

export FullDynamicState, DynamicState

"""
Stores the relevant Dynamic variables of the structure during the displacements iteration.
"""
struct FullDynamicState{DU<:AbstractVector,U<:AbstractVector,
                        FE<:AbstractVector,FI<:AbstractVector,K<:AbstractMatrix,
                        E<:Dictionary,S<:Dictionary} <: AbstractDynamicState
    "Free degrees of freedom."
    free_dofs::Vector{Dof}
    "Displacements vector increment."
    ΔUᵏ::DU
    "Displacements vector."
    Uᵏ::U
    "Velocity vector."
    Udotᵏ::U
    "Acceleration vector."
    Udotdotᵏ::U
    "External forces vector."
    Fₑₓₜᵏ::FE
    "Internal forces vector."
    Fᵢₙₜᵏ::FI
    "Inertial forces vector."
    Fᵢₙₑᵏ::FI
    "Damping forces vector."
    Fᵥᵢₛᵏ::FI
    "Stiffness tangent matrix."
    Kᵏ::K
    "Mass Matrix tangent matrix."
    Mᵏ::K
    "Damping Matrix tangent matrix."
    Cᵏ::K
    "System tangent matrix."
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
    "Linear system cache"
    linear_system::LinearSolve.LinearCache
    function FullDynamicState(fdofs::Vector{Dof},
                              ΔUᵏ::DU, Uᵏ::U, Udotᵏ::U, Udotdotᵏ::U
                              Fₑₓₜᵏ::FE, Fᵢₙₜᵏ::FI,Fᵢₙₑᵏ::FI, Fᵥᵢₛᵏ::FI,
                              Kᵏ::K,Mᵏ::K, Cᵏ::K, Kₛᵏ::K, res_forces::DU,
                              ϵᵏ::E, σᵏ::S,
                              assembler::Assembler,
                              iter_state::ResidualsIterationStep,
                              linear_system::LinearSolve.LinearCache) where {DU,U,FE,FI,K,E,S}
        # Check dimensions
        @assert length(ΔUᵏ) == length(fdofs) == length(res_forces)
        @assert begin
            size(Kₛᵏ, 1) == size(Kₛᵏ, 2) == size(Mᵏ,1) == size(Mᵏ,2) == size(Cᵏ,1) == size(Cᵏ,2) ==  size(Kᵏ,1) == size(Kᵏ,2)  == length(Fᵢₙₜᵏ) == length(Fₑₓₜᵏ) == length(Uᵏ)
        end
        # Initialize linear system K.ΔU = R
        new{DU,U,FE,FI,K,E,S}(
            fdofs,
        ΔUᵏ, Uᵏ, Udotᵏ, Udotdotᵏ,
        Fₑₓₜᵏ, Fᵢₙₜᵏ,Fᵢₙₑᵏ, Fᵥᵢₛᵏ,
        Kᵏ,Mᵏ, Cᵏ, Kₛᵏ, res_forces, ϵᵏ, σᵏ,
        assembler, iter_state, linear_system)
    end
end

"Default constructor for Dynamic state given an structure and iteration state."
function FullDynamicState(s::AbstractStructure,
                          iter_state::ResidualsIterationStep=ResidualsIterationStep())
    n_dofs = num_dofs(s)
    n_fdofs = num_free_dofs(s)

    Uᵏ = zeros(n_dofs)
    Udotᵏ = zeros(n_dofs)
    Udotdotᵏ = zeros(n_dofs)
    ΔUᵏ = zeros(n_fdofs)

    Fₑₓₜᵏ = zeros(n_dofs)
    Fᵢₙₜᵏ = zeros(n_dofs)
    Fᵢₙₑᵏ = zeros(n_dofs)
    Fᵥᵢₛᵏ = zeros(n_dofs)


    Kᵏ = spzeros(n_dofs, n_dofs)
    Mᵏ = spzeros(n_dofs, n_dofs)
    Cᵏ = spzeros(n_dofs, n_dofs)


    Kₛᵏ = spzeros(n_dofs, n_dofs)
    res_forces = zeros(n_fdofs)

    # Initialize pairs strains
    ϵᵏ = dictionary([Pair(e, Symmetric(Matrix{Float64}(undef, (3, 3)))) for e in elements(s)])
    σᵏ = dictionary([Pair(e, Matrix{Float64}(undef, (3, 3))) for e in elements(s)])

    cache = dictionary(nameof(T) => elements_cache(T) for T in subtypes(AbstractElement))
    assemblerᵏ = Assembler(s, cache)

    fdofs = free_dofs(s)
    linear_system = init(LinearProblem(Kₛᵏ[fdofs, fdofs], res_forces))

    FullDynamicState(fdofs, ΔUᵏ, Uᵏ,Udotᵏ,Udotdotᵏ Fₑₓₜᵏ, Fᵢₙₜᵏ, Kₛᵏ, res_forces, ϵᵏ, σᵏ, assemblerᵏ, iter_state,
                     linear_system)
end

̈U = 2
function Base.show(io::IO, sc::FullDynamicState)
    nu = length(sc.Uᵏ)
    K = sc.Kₛᵏ
    s = size(K)
    println("• FullDynamicState with $nu-dofs displacements vector Uᵏ " *
            "and $(s[1]) × $(s[2]) tangent matrix Kₛᵏ with $(length(K.nzval)) stored entries.")
end

"Update and return the current residual forces of the Dynamic state."
function residual_forces!(sc::FullDynamicState)
    return sc.res_forces .= view(external_forces(sc), free_dofs(sc)) -
                            view(internal_forces(sc), free_dofs(sc))
end

"Return the current system tangent matrix form the Dynamic state ."
tangent_matrix(sc::FullDynamicState) = sc.Kₛᵏ

"Reset the Dynamic state assembled magnitudes and the iteration state."
function reset!(state::FullDynamicState)
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

struct DynamicState{U<:AbstractVector,E<:Dictionary,S<:Dictionary} <: AbstractDynamicState
    "Displacements vector."
    Uᵏ::U
    "Vector with straings for each element."
    ϵᵏ::E
    "Vector with stresses for each element."
    σᵏ::S
    function DynamicState(Uᵏ::U, ϵᵏ::E, σᵏ::S) where {U,E,S}
        new{U,E,S}(Uᵏ, ϵᵏ, σᵏ)
    end
end

function Base.show(io::IO, sc::DynamicState)
    nu = length(sc.Uᵏ)
    println("• DynamicState with $nu-dofs displacements vector Uᵏ.")
end

end # module
