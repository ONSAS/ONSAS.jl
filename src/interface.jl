
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


mutable struct CrossSection
    type::String
    diameter
    width_y
    width_z
end

function CrossSection(type, params)
    if cmp( type, "square") == 0
        return CrossSection( type, nothing, params[1], params[1] )
    elseif cmp( type, "circle") == 0
        return CrossSection( type, params[1], nothing, nothing )
    end
end


"""
Struct with information of the geometry of the element:
 - element_type: a string with: `node`, `truss`, `frame``, `triangle` or `tetrahedron`.
"""
struct ElementGeometry
    element_type::String
    cross_section::CrossSection
end

mutable struct BoundaryCondsData
    DirichletNodalDOFs::Vector{Int}
    DirichletNodalVals::Vector{Float64}
    NeumannNodalDOFs::Vector{Int}
    NeumannNodalVals::Vector{Float64}
end

struct BCsIndexes
    diriDofs::Vector{Int}
    neumDofs::Vector{Int}
end

struct SystemMatrix{T}
    matrix::T
end
