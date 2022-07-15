
# ======================================================================
# material
# ======================================================================

"""
Struct with information of the material.
"""
mutable struct Material{T}
    type::String
    constitutive_params::Vector{T}
    density::Float64 
end

function Material(type, constitutive_params)
    return Material(type, constitutive_params,0.0)
end






# ======================================================================
# geometry
# ======================================================================

mutable struct CrossSection
    type::String
    diameter
    width_y
    width_z
    area::Float64
    inertia_x::Float64
    inertia_y::Float64
    inertia_z::Float64
end

# constructor with missing fields
function CrossSection(type; width_y = nothing, width_z = nothing )
    if cmp( type, "square") == 0 || cmp( type, "rectangle") == 0
        
        # set aditional width for square section case
        if cmp( type, "square") == 0
            if isnothing( width_y)
                width_y = width_z
            else
                width_z = width_y
            end
        end

        area      = width_y * width_z 

        inertia_y = width_y * width_z^3 / 12
        inertia_z = width_z * width_y^3 / 12

        # torsional constant from table 10.1 from Roark's Formulas for Stress and Strain 7th ed.
        a = 0.5 * max( width_y, width_z )
        b = 0.5 * min( width_y, width_z )
            
        inertia_x = a * b^3 * ( 16/3. - 3.36 * b/a * ( 1. - b^4 / ( 12*a^4 ) ) )
        
        return CrossSection( type, nothing, width_y, width_z, area, inertia_x, inertia_y, inertia_z )
    else
      error(" cross section type ", type, " not implemented yet, please create an issue!")
    end
end


"""
Struct with information of the geometry of the element:
 - element_type: a string with: `node`, `truss`, `frame``, `triangle` or `tetrahedron`.
"""
struct Geometry
    type::String
    cross_section
end

# constructor with missing fields
function Geometry( type::String )
    return Geometry( type, nothing )
end


# ======================================================================
# Boundary Conditions
# ======================================================================

mutable struct BoundaryCondition
    imposed_disp_dofs::Vector{Int}
    imposed_disp_vals::Vector{Float64}
    user_load_function::String
end

# ======================================================================
# Initial Conditions
# ======================================================================

struct InitialCondition
    dofs::Vector{Int}
    vals::Vector{Float64}
end


# ======================================================================
# Mesh
# ======================================================================

struct Mesh
    nodal_coords::Matrix{Float64}
    elem_nodal_connec::Vector{Vector{Int}}
    MGBI_mat::Matrix{Int}
    MGBI_vec::Vector{Int}
end


# ======================================================================
# AnalysisSettings
# ======================================================================

struct AnalysisSettings

    method::String
    delta_time::Float64
    final_time::Float64

#    delta_time > final_time && error("delta_time must be lower than final_time")

    stop_tol_disps::Float64
    stop_tol_force::Float64
    stop_tol_iters::Int
end

function AnalysisSettings( method::String, delta_time::Float64, final_time::Float64 )
    return AnalysisSettings( method, delta_time, final_time, 1.0e-6, 1.0e-6, 20 )
end

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
displacements( sol::ModelSolution ) = sol.U

function unwrap( sol::ModelSolution )
    return sol.time, sol.U, sol.Udot, sol.Udotdot 
end


struct ModelProperties
    materials::Vector{Material}
    geometries::Vector{Geometry}
    boundary_conditions::Vector{BoundaryCondition}
    neum_dofs::Vector{Int}
    mesh::Mesh
    analysis_settings::AnalysisSettings
end




struct SystemMatrix{T}
    matrix::T
end






