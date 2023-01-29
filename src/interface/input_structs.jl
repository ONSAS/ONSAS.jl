
include("./Materials.jl")
@reexport using .Materials


include("./BoundaryConditions.jl")
@reexport using .BoundaryConditions

include("./CrossSections.jl")
@reexport using .CrossSections

include("./Meshes.jl")
@reexport using .Meshes


"""
    AbstractElement

Abstract type to define finite element types.\\

Available types:

* `Node`: 
* `Truss`: 
* `Frame`: 
* `Triangle`: 
* `Tetrahedron`: 


"""
abstract type AbstractElement end

# struct Node <: AbstractElement end
# struct Truss <: AbstractElement end
struct Frame <: AbstractElement end
struct Triangle <: AbstractElement end
struct Tetrahedron <: AbstractElement end

"""
Struct with information of the geometry of the element:
 - element_type: a struct with: `node`, `truss`, `frame``, `triangle` or `tetrahedron`.
"""
struct Geometry
    type::AbstractElement
    cross_section
    # define constructor with no cross section by default
    function Geometry(type::AbstractElement, cross_section=nothing)
        new(type, CrossSection())
    end
end



# ======================================================================
# Dofs Boundary Conditions - Springs and imposed displacements (zero & nonzero)
# ======================================================================
abstract type AbstractDofs end


# ======================================================================
# Initial Conditions
# ======================================================================
struct InitialCondition
    dofs::Vector{Integer}
    vals::Vector{Float64}
end


# ======================================================================
# Mesh
# ======================================================================
# struct Mesh
#     nodal_coords::Matrix{Float64}
#     elem_nodal_connec::Vector{Vector{Int64}}
#     MGBI_mat::Matrix{Int64}
#     MGBI_vec::Vector{Int64}
# end


# ======================================================================
# AnalysisSettings
# ======================================================================
"""
    AbstractAlgorithm

Abstract type to define numerical method algorithms.

"""
abstract type AbstractAlgorithm end

"""
    AnalysisSettings

Struct to define convergence tolerances.

"""
struct ConvergenceSettings

    # method::String
    # delta_time::Float64
    # final_time::Float64

    #delta_time > final_time && error("delta_time must be lower than final_time")

    stop_tol_disps::Float64
    stop_tol_force::Float64
    stop_tol_iters::Int

    # function AnalysisSettings(method::String, delta_time::Float64, final_time::Float64,
    #     stop_tol_disps=1e-6::Float64, stop_tol_force=1e-6::Float64, stop_tol_iters=20::Integer)
    #     new(method, delta_time, final_time, stop_tol_disps, stop_tol_force, stop_tol_iters)
    # end

    function ConvergenceSettings(stop_tol_disps=1e-6::Float64, stop_tol_force=1e-6::Float64, stop_tol_iters=20::Int)
        new(stop_tol_disps, stop_tol_force, stop_tol_iters)
    end

end
