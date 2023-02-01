"""
Module defining structure entities interface.
"""
module StructuralModel

using ..Utils: set_label!
using Reexport: @reexport

export StructuralMaterials, StructuralElements

# ======================
# Model definition
# ======================

""" Structural materials.
A `StructuralMaterials` is a collection of `Materials` defining the material models of the structure.
### Fields:
- `materials` -- Stores material models` 
"""
struct StructuralMaterials{M}
    materials::Vector{M}
end

function StructuralMaterials(mats_dict::Dict{String,M}) where {M}

    vmats = Vector{M}(undef, length(mats_dict))
    for (l, m) in mats_dict
        set_label!(m, l)
        push!(vmats, m)
    end

    return StructuralMaterials(vmats)
end

""" Structural elements.
A `StructuralElements` is a collection of `Elements` defining the elements types of the structure.
### Fields:
- `elements` -- Stores element types. 
"""
struct StructuralElements{E}
    materials::Vector{E}
end

function StructuralElements(elems_dict::Dict{String,E}) where {E}

    velems = Vector{E}(undef, length(elems_dict))
    for (l, e) in elems_dict
        set_label!(e, l)
        push!(velems, e)
    end

    return StructuralElements(velems)
end


#=
""" Structural boundary conditions.
A `StructuralBoundaryConditions` is a collection of `BoundaryConditions` defining the boundary conditions of the structure.
### Fields:
- `displacements_bc` -- Stores displacement boundary conditions. 
- `displacements_bc` -- Stores loads boundary conditions. 
"""
struct StructuralBoundaryConditions{E}
    materials::Vector{E}
end

function StructuralElements(elems_dict::Dict{String,E}) where {E}

    velems = Vector{E}(undef, length(elems_dict))
    for (l, e) in elems_dict
        set_label!(e, l)
        push!(velems, e)
    end

    return StructuralElements(velems)
end

=#
end # module
