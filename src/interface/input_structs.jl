
# ======================================================================
# material
# ======================================================================

"""
Struct with information of the material.
"""
mutable struct Material{T}
    type::String
    constitutive_params::Vector{T}
    density::Union{Float64, Nothing}
    # define constructor with  no density by default
    function Material(type::String, constitutive_params::Vector{T}, density=nothing::Union{Float64, Nothing}) where T
        new{T}(type, constitutive_params, density) # default density: zero
    end
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
    # define constructor with no cross section by default
    function Geometry( type::String, cross_section = nothing  )
        new( type, cross_section )
    end
end

# ======================================================================
# boundary conditions
# ======================================================================
mutable struct BoundaryCondition
    imposed_disp_dofs::Vector{Int}
    imposed_disp_vals::Vector{Float64}
    user_load_function
end
# constructor with missing fields
function BoundaryCondition( imposed_disp_dofs, imposed_disp_vals )
    return BoundaryCondition( imposed_disp_dofs, imposed_disp_vals, nothing )
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
    elem_nodal_connec::Vector{Vector{Int64}}
    MGBI_mat::Matrix{Int64}
    MGBI_vec::Vector{Int64}
end


# ======================================================================
# AnalysisSettings
# ======================================================================

struct AnalysisSettings

    method::String
    delta_time::Float64
    final_time::Float64

    #delta_time > final_time && error("delta_time must be lower than final_time")

    stop_tol_disps::Float64
    stop_tol_force::Float64
    stop_tol_iters::Int

    function AnalysisSettings(  method::String, delta_time::Float64, final_time::Float64, 
                                stop_tol_disps=1e-6::Float64, stop_tol_force=1e-6::Float64, stop_tol_iters=20::Integer )
        new( method, delta_time, final_time, stop_tol_disps, stop_tol_force, stop_tol_iters )
    end

end
