include("StructuralModel.jl")
@reexport using .StructuralModel

include("StructuralAnalyses.jl")
@reexport using .StructuralAnalyses




# # ======================================================================
# # Mesh
# # ======================================================================
# # struct Mesh
# #     nodal_coords::Matrix{Float64}
# #     elem_nodal_connec::Vector{Vector{Int64}}
# #     MGBI_mat::Matrix{Int64}
# #     MGBI_vec::Vector{Int64}
# # end


# # ======================================================================
# # AnalysisSettings
# # ======================================================================
# """
#     AbstractAlgorithm

# Abstract type to define numerical method algorithms.

# """
# abstract type AbstractAlgorithm end

# """
#     AnalysisSettings

# Struct to define convergence tolerances.

# """
# struct ConvergenceSettings

#     # method::String
#     # delta_time::Float64
#     # final_time::Float64

#     #delta_time > final_time && error("delta_time must be lower than final_time")

#     stop_tol_disps::Float64
#     stop_tol_force::Float64
#     stop_tol_iters::Int

#     # function AnalysisSettings(method::String, delta_time::Float64, final_time::Float64,
#     #     stop_tol_disps=1e-6::Float64, stop_tol_force=1e-6::Float64, stop_tol_iters=20::Integer)
#     #     new(method, delta_time, final_time, stop_tol_disps, stop_tol_force, stop_tol_iters)
#     # end

#     function ConvergenceSettings(stop_tol_disps=1e-6::Float64, stop_tol_force=1e-6::Float64, stop_tol_iters=20::Int)
#         new(stop_tol_disps, stop_tol_force, stop_tol_iters)
#     end

# end






# mutable struct ModelSolution
#     time::Float64
#     U::Vector{Float64}       # displacements
#     Udot::Vector{Float64}    # velocities
#     Udotdot::Vector{Float64} # accelerations
#     system_matrix
#     system_rhs
# end

# # function ModelSolution( time, U, Udot, Udotdot )
# #     return ModelSolution(time, U, Udot, Udotdot, nothing, nothing )
# # end

# """
#     displacements(sol::ModelSolution)
# Return the vector of displacements of the given solution.
# """
# displacements(sol::ModelSolution) = sol.U

# function unwrap(sol::ModelSolution)
#     return sol.time, sol.U, sol.Udot, sol.Udotdot
# end


# # struct ModelProperties
# #     Materials::Vector{<:AbstractMaterial}
# #     Geometries::Vector{Geometry}
# #     LoadsBC::Vector{LoadsBoundaryCondition}
# #     DofsBC::Vector{DispsBoundaryCondition}
# #     Mesh::Mesh
# #     ConvSettings::ConvergenceSettings
# #     Algorithm::AbstractAlgorithm
# #     neum_dofs::Vector{Int}
# # end
