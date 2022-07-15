
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
    Area::Float64
    Inertia_x::Float64
    Inertia_y::Float64
    Inertia_z::Float64
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

        Area      = width_y * width_z 

        Inertia_y = width_y * width_z^3 / 12
        Inertia_z = width_z * width_y^3 / 12

        # torsional constant from table 10.1 from Roark's Formulas for Stress and Strain 7th ed.
        a = 0.5 * max( width_y, width_z )
        b = 0.5 * min( width_y, width_z )
            
        Inertia_x = a * b^3 * ( 16/3. - 3.36 * b/a * ( 1. - b^4 / ( 12*a^4 ) ) )
        
        return CrossSection( type, nothing, width_y, width_z, Area, Inertia_x, Inertia_y, Inertia_z )
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


struct AnalysisSettings
    numerical_method::String
    delta_time::Float64
    final_time::Float64
end

struct ModelSolution
    time::Float64
    displacements::Vector{Float64}
    velocities::Vector{Float64}
    accelerations::Vector{Float64}
end


struct ModelProperties
    mesh::Mesh
    analysis_settings::AnalysisSettings
end




struct SystemMatrix{T}
    matrix::T
end






