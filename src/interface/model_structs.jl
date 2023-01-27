

mutable struct ModelSolution
    time::Float64
    U::Vector{Float64}       # displacements
    Udot::Vector{Float64}    # velocities
    Udotdot::Vector{Float64} # accelerations
    system_matrix
    system_rhs
end

# function ModelSolution( time, U, Udot, Udotdot )
#     return ModelSolution(time, U, Udot, Udotdot, nothing, nothing )
# end

"""
    displacements(sol::ModelSolution)
Return the vector of displacements of the given solution.
"""
displacements(sol::ModelSolution) = sol.U

function unwrap(sol::ModelSolution)
    return sol.time, sol.U, sol.Udot, sol.Udotdot
end


# struct ModelProperties
#     Materials::Vector{<:AbstractMaterial}
#     Geometries::Vector{Geometry}
#     LoadsBC::Vector{LoadsBoundaryCondition}
#     DofsBC::Vector{DispsBoundaryCondition}
#     Mesh::Mesh
#     ConvSettings::ConvergenceSettings
#     Algorithm::AbstractAlgorithm
#     neum_dofs::Vector{Int}
# end




struct SystemMatrix{T}
    matrix::T
end
