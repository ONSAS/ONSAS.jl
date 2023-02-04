"""
Module defining structure entities interface.
"""
module StructuralModel

using Reexport: @reexport
using StaticArrays: @MVector
using SparseArrays: SparseMatrixCSC

using ..Materials: AbstractMaterial
using ..CrossSections: AbstractCrossSection
using ..Elements
using ..BoundaryConditions: AbstractBoundaryCondition, AbstractLoadBoundaryCondition, AbstractDisplacementBoundaryCondition
using ..Utils: label, set_label!
using ..Meshes: AbstractMesh, element_nodes

@reexport import ..Meshes: elements
@reexport import ..Elements: nodes
@reexport import ..Utils: external_forces, internal_forces, internal_tangents, displacements

export AbstractStructuralMEBI, sets, add_set!
export StructuralMaterials, materials
export StructuralElements
export StructuralBoundaryConditions, disp_bcs, load_bcs
export StructuralInitialConditions
export AbstractStructure, Structure, mesh, current_state
export AbstractStructuralState, StaticState

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
    new_set::Dict{String,<:Set{<:AbstractIndex}}
)
    merge!(sets(mebi), new_set)
end

#TODO: Check not repeated keys 

""" Structural materials.
A `StructuralMaterials` is a collection of `Materials` defining the material models of the structure.
### Fields:
- `vec_mats` -- Stores material models`.
- `sets`     -- Maps a material `String` (corresponding to a material type) into element ids. 
"""
struct StructuralMaterials{M} <: AbstractStructuralMEBI
    vec_mats::Vector{M}
    sets::Dict{String,<:Set}
    function StructuralMaterials(
        mats_dict::Dict{String,M},
        sets::Dict{String,<:Set}=Dict{String,Set}()
    ) where {M<:AbstractMaterial}

        vmats = Vector{M}(undef, 0)
        for (l, m) in mats_dict
            set_label!(m, l)
            push!(vmats, m)
        end

        new{M}(vmats, sets)
    end
end

materials(sm::StructuralMaterials) = sm.vec_mats

function Base.getindex(sm::StructuralMaterials, mat_label::L) where {L<:Union{String,Symbol}}
    mat_vec = materials(sm)
    mat = filter(m -> label(m) == Symbol(mat_label), mat_vec)
    isempty(mat) && @warn("The label $mat_label was not found among the structural materials")

    return first(mat)
end

"Returns material type for a given index"
function _find_material(sm::StructuralMaterials, i_e::ElementIndex)
    # Find element material
    for (label_mat_set, mat_set) in sets(sm)
        if i_e ∈ mat_set
            mat_elem = sm[label_mat_set]
            return mat_elem
        end
    end
    nothing
end

""" Structural elements.
A `StructuralElements` is a collection of `Elements` defining the elements types of the structure.
### Fields:
- `vec_elems` -- Stores element types. 
- `sets`     -- Maps an element `String` (corresponding to a element type) into element ids. 
"""
struct StructuralElements{E} <: AbstractStructuralMEBI
    vec_elems::Vector{E}
    sets::Dict{String,<:Set}
    function StructuralElements(
        elems_dict::Dict{String,E},
        sets::Dict{String,<:Set}=Dict{String,Set{AbstractIndex}}()
    ) where {E<:AbstractElement}

        velems = Vector{E}(undef, 0)
        for (l, e) in elems_dict
            set_label!(e, l)
            push!(velems, e)
        end

        new{E}(velems, sets)
    end
end

elements(se::StructuralElements) = se.vec_elems

function Base.getindex(se::StructuralElements, elem_label::L) where {L<:Union{String,Symbol}}
    elem_vec = elements(se)
    elems = filter(e -> label(e) == Symbol(elem_label), elem_vec)
    isempty(elems) && @warn "The label $elem_label was not found among the structural elements"
    return first(elems)
end

"Returns material type for a given index"
function _find_element_type(se::StructuralElements, i::ElementIndex)
    # Find element elementype
    for (label_elem_set, elem_set) in sets(se)
        if i ∈ elem_set
            elemtype = se[label_elem_set]
            return elemtype
        end
    end
    nothing
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
    sets::Dict{String,<:Set}
    function StructuralBoundaryConditions(
        bc_dict::Dict{String,B},
        sets::Dict{String,<:Set}=Dict{String,Set{AbstractIndex}}()
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

disp_bcs(se::StructuralBoundaryConditions) = se.displacements_bc
load_bcs(se::StructuralBoundaryConditions) = se.loads_bc

function Base.getindex(bcs::StructuralBoundaryConditions, bc_label::L) where {L<:Union{String,Symbol}}
    disp_bcs_vec = disp_bcs(bcs)
    disp_bcs_label_vec = filter(dbc -> label(dbc) == Symbol(bc_label), disp_bcs_vec)
    isempty(disp_bcs_label_vec) || return first(disp_bcs_label_vec)
    load_bcs_vec = load_bcs(bcs)
    load_bcs_label_vec = filter(lbc -> label(lbc) == Symbol(bc_label), load_bcs_vec)
    isempty(load_bcs_label_vec) || return first(load_bcs_label_vec)
    return @warn "The label $bc_label was not found among the boundary conditions"
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
    sets::Dict{String,Set{AbstractIndex}}
    function StructuralInitialConditions(
        ic_dict::Dict{String,I},
        sets::Dict{String,Set{AbstractIndex}}=Dict{String,Set{AbstractIndex}}()
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
# Structural state
# ======================


""" Abstract supertype to define a new structural state.
structure(s::AbstractStructuralState) = s.s
**Common methods:**
* [`displacements`](@ref)
* [`internal_forces`](@ref)
* [`internal_tangents`](@ref)
* [`external_forces`](@ref)
"""
abstract type AbstractStructuralState end

displacements(sc::AbstractStructuralState) = sc.Uᵏ
internal_forces(sc::AbstractStructuralState) = sc.Fᵢₙₜᵏ
internal_tangents(sc::AbstractStructuralState) = sc.Kₛᵏ
external_forces(sc::AbstractStructuralState) = sc.Fₑₓₜᵏ

"""
An `StaticState` object facilitates the process of storing the relevant static variables of the structure. 
### Fields:
- `Uᵏ`    -- stores displacements vector.
- `Fₑₓₜᵏ` -- stores external forces vector.
- `Fᵢₙₜᵏ` -- stores internal forces vector.
- `Kₛᵏ`   -- stiffness tangent matrix of the structure. 
"""
struct StaticState <: AbstractStructuralState
    Uᵏ::AbstractVector
    Fₑₓₜᵏ::AbstractVector
    Fᵢₙₜᵏ::AbstractVector
    Kₛᵏ::AbstractMatrix
end

#TODO: Add tangent matrix of the external forces vector

"Returns a default static case for a given mesh."
function StaticCase(m::AbstractMesh)
    n_dofs = length(dofs(m))
    Uᵏ = @MVector zeros(n_dofs)
    Fₑₓₜᵏ = similar(Uᵏ)
    Fᵢₙₜᵏ = similar(Uᵏ)
    Kₛᵏ = SparseMatrixCSC(zeros(n_dofs, n_dofs))
    StaticState(Uᵏ, Fₑₓₜᵏ, Fᵢₙₜᵏ, Kₛᵏ)
end

# ==========
# Structure
# ==========

""" Abstract supertype to define a new structure.
An `AbstractStructure` object facilitates the process of assigning materials, 
elements, initial and boundary conditions to the mesh.
**Common methods:**
* [`nodes`](@ref)
* [`elements`](@ref)
* [`mesh`](@ref)
"""
abstract type AbstractStructure{D,M,E,B,I} end

"""
An `Structure` object facilitates the process of assembling and creating the structural analysis. 
### Fields:
- `mesh`        -- Stores the structure mesh. 
- `materials`   -- Stores the types of material models considered in the structure. 
- `elements`    -- Stores the types of elements considered in the structure.
- `bcs`         -- Stores the types of boundary conditions in the structure.
- `ics`         -- Stores the types of initial conditions in the structure.
- `ics`         -- Stores the current state of the structure.
"""
mutable struct Structure{D,M,E,B,I} <: AbstractStructure{D,M,E,B,I}
    mesh::AbstractMesh{D}
    materials::StructuralMaterials{M}
    elements::StructuralElements{E}
    bcs::StructuralBoundaryConditions{B}
    ics::StructuralInitialConditions{I}
    state::AbstractStructuralState
    function Structure(
        mesh::AbstractMesh{D},
        mats::StructuralMaterials{M},
        elems::StructuralElements{E},
        bcs::StructuralBoundaryConditions{B},
        init::StructuralInitialConditions{I}=StructuralInitialConditions{Nothing}()
    ) where {D,M,E,B,I}


        # Create elements and push them into the mesh
        _create_elements!(mesh, mats, elems)
        # Apply bc
        _apply_bcs!(mesh, bcs)

        # Default structural state
        s_case = StaticCase(mesh)

        return new{D,M,E,B,I}(mesh, mats, elems, bcs, init, s_case)
    end
end

"Creates a new element given a pre element (without material) defined in the input"
function _create_elements!(
    mesh::AbstractMesh,
    mats::StructuralMaterials,
    elems::StructuralElements
)

    for (i, e_nodes) in enumerate(element_nodes(mesh)) # loop over elements

        e_mat = _find_material(mats, ElementIndex(i))
        pre_element = _find_element_type(elems, ElementIndex(i))
        e = _create_element(pre_element, e_mat, e_nodes)
        push!(mesh, e)
    end

end

"Function to assign material to a truss with a given pre_element"
function _create_element(pre_element::Truss, m::AbstractMaterial, nodes::AbstractVector{<:AbstractNode})
    element_type(pre_element)(nodes, m, geometry(pre_element), pre_element.label)
end

"Applies nodal and element boundary conditions to the mesh"
function _apply_bcs!(mesh::AbstractMesh, bcs::StructuralBoundaryConditions)

    for (label, set) in sets(bcs)
        bc = bcs[label]
        for index in set
            push!(mesh[index], bc)
        end
    end

end

"Returns the mesh"
mesh(s::Structure) = s.mesh

nodes(s::Structure) = nodes(mesh(s))

elements(s::Structure) = elements(mesh(s))

Base.getindex(s::Structure, i::AbstractIndex) = mesh(s)[i]

"Returns the current structural state"
current_state(s::Structure) = s.state



end # module
