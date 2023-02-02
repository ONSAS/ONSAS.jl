"""
Module defining structure entities interface.
"""
module StructuralModel

using ..Materials: AbstractMaterial
using ..Elements: AbstractElement
using ..BoundaryConditions: AbstractBoundaryCondition, AbstractLoadBoundaryCondition, AbstractDisplacementBoundaryCondition
using ..Utils: label, set_label!
using ..Meshes: AbstractMesh, element_nodes
using Reexport: @reexport

export StructuralDefinitionStruct, StructuralMaterials, StructuralElements,
    StructuralBoundaryConditions, StructuralInitialConditions, Structure

export sets, add_set!

# ======================
# Model definition
# ======================

""" Abstract supertype for all structural definitions.

An `AbstractStructuralMEBI` object facilitates the process of defining sets of materials, elements, initial and boundary conditions.

**Common methods:**

* [`sets`](@ref)
* [`add_set!`](@ref)
"""
abstract type AbstractStructuralMEBI end

"Returns sets defined."
sets(mebi::AbstractStructuralMEBI) = mebi.sets

"Adds a new set to the sets."
function add_set!(
    mebi::AbstractStructuralMEBI,
    new_set::Dict{String,Set{Int}}
)
    merge!(sets(mebi), new_set)
end

#Add getindex mat["mat1"] = material 1

""" Structural materials.
A `StructuralMaterials` is a collection of `Materials` defining the material models of the structure.
### Fields:
- `vec_mats` -- Stores material models`.
- `sets`     -- Maps a material `String` (corresponding to a material type) into element ids. 
"""
struct StructuralMaterials{M} <: AbstractStructuralMEBI
    vec_mats::Vector{M}
    sets::Dict{String,Set{Int}}
    function StructuralMaterials(
        mats_dict::Dict{String,M},
        sets::Dict{String,Set{Int}}=Dict{String,Set{Int}}()
    ) where {M<:AbstractMaterial}

        vmats = Vector{M}(undef, 0)
        for (l, m) in mats_dict
            set_label!(m, l)
            push!(vmats, m)
        end

        new{M}(vmats, sets)
    end
end

""" Structural elements.
A `StructuralElements` is a collection of `Elements` defining the elements types of the structure.
### Fields:
- `vec_elems` -- Stores element types. 
- `sets`     -- Maps an element `String` (corresponding to a element type) into element ids. 
"""
struct StructuralElements{E} <: AbstractStructuralMEBI
    vec_elems::Vector{E}
    sets::Dict{String,Set{Int}}
    function StructuralElements(
        elems_dict::Dict{String,E},
        sets::Dict{String,Set{Int}}=Dict{String,Set{Int}}()
    ) where {E<:AbstractElement}

        velems = Vector{E}(undef, 0)
        for (l, e) in elems_dict
            set_label!(e, l)
            push!(velems, e)
        end

        new{E}(velems, sets)
    end
end


""" Structural boundary conditions.
A `StructuralBoundaryConditions` is a collection of `BoundaryConditions` defining the boundary conditions of the structure.
### Fields:
- `displacements_bc` -- Stores displacement boundary conditions. 
- `loads_bc` -- Stores loads boundary conditions. 
- `loads_bc` -- Stores loads boundary conditions. 
- `sets`     -- Maps a boundary condition `String` (corresponding to a bc type) into element ids. 
"""
struct StructuralBoundaryConditions{B} <: AbstractStructuralMEBI
    displacements_bc::Vector{B}
    loads_bc::Vector{B}
    sets::Dict{String,Set{Int}}
    function StructuralBoundaryConditions(
        bc_dict::Dict{String,B},
        sets::Dict{String,Set{Int}}=Dict{String,Set{Int}}()
    ) where {B<:AbstractBoundaryCondition}

        v_dbc = Vector{B}(undef, 0)
        v_lbc = Vector{B}(undef, 0)
        for (l, b) in bc_dict
            set_label!(b, l)
            if b isa AbstractLoadBoundaryCondition
                push!(v_lbc, b)
            elseif b isa AbstractDisplacementBoundaryCondition
                push!(v_dbc, b)
            end
        end

        new{B}(v_dbc, v_lbc, sets)
    end
end

# TODO: Implement initial conditions struct
struct StructuralInitialConditions{I} <: AbstractStructuralMEBI end

#=

const D = DisplacementInitialCondition
const V = VelocityInitialCondition
const A = AccelerationInitialCondition


""" Structural initial conditions.
A `StructuralInitialConditions` is a collection of `InitialConditions` defining the initial conditions of the structure.
### Fields:
- `displacements`-- Stores displacement boundary conditions. 
- `velocity`     -- Stores displacement boundary conditions. 
- `acceleration` -- Stores displacement boundary conditions. 
- `sets`         -- Maps an initial condition `String` (corresponding to a ic type) into elements ids. 
"""
struct StructuralInitialConditions{D,V,A} <: AbstractStructuralMEBI
    displacements::Vector{D}
    velocity::Vector{V}
    acceleration::Vector{A}
    sets::Dict{String,Set{Int}}
    function StructuralInitialConditions(
        ic_dict::Dict{String,I},
        sets::Dict{String,Set{Int}}=Dict{String,Set{Int}}()
    ) where {I<:AbstractElement}

        v_dic = Vector{D}(undef, 0)
        v_vic = Vector{V}(undef, 0)
        v_aic = Vector{A}(undef, 0)
        for (l, i) in ic_dict
            set_label!(i, l)
            if b isa D
                push!(v_dic, i)
            elseif b isa V
                push!(v_vic, i)
            elseif b isa A
                push!(v_aic, i)
            end
        end

        new{D,V,A}(v_dic, v_vic, v_aic, sets)
    end
end
=#


# ======================
# Structure
# ======================


""" Abstract supertype to define a new structure.

An `AbstractStructure` object facilitates the process of assigning materials, 
elements, initial and boundary conditions to the mesh.

**Common methods:**

* [`nodes`](@ref)
* [`elements`](@ref)

* [`materials`](@ref)
* [`boundary_conditions`](@ref)
* [`initial_conditions`](@ref)

* [`displacements`](@ref)
* [`velocities`](@ref)
* [`accelerations`](@ref)

* [`tangent_matrices`](@ref)
* [`forces_vectors`](@ref)

"""
abstract type AbstractStructure{M,E,B,I} end


"""
An `Structure` object facilitates the process of asse 
### Fields:
- `mats` -- Stores a list with each material defined object. 
- `elements` -- Stores a list with elements. 
- `bcs` -- Stores a list with different boundary conditions. 
- `initial` -- Stores a list with initial conditions. 
"""
struct Structure{D,M,E,B,I}
    mats::AbstractVector{M}
    elements::AbstractVector{E}
    bcs::AbstractVector{B}
    initial::AbstractVector{I}
    mesh::AbstractMesh{D}
end


function Structure(
    mesh::AbstractMesh{D},
    mats::StructuralMaterials{M},
    elements::StructuralElements{E},
    bcs::StructuralBoundaryConditions{B},
    init::StructuralInitialConditions{I}=StructuralInitialConditions{Nothing}()
) where {D,M,E,B,I}

    Main.@infiltrate

    # set_mats = set(mats)
    # set_elements = set(elements)
    # set_elements = set(bcs)
    for (e_nodes, i) in enumerate(element_nodes(mesh)) # loop over elements
        for (label_mat_set, set) in mats
            if i ∈ set
                break
            end
            mat_label = label_mat_set
            mat_e = set(materials)[mat_label]
            # element_material_label = [mat  if i ∈ set]            
        end

    end
    # Structure{D,M,E,B,I}(mats)
    return nothing
end



end # module
