"Module defining the structural boundary conditions assigned to each entitiy. This data type
includes:

- `AbstractBoundaryCondition`s assigned to `Node`s
- `AbstractBoundaryCondition`s assigned to `Face`s
- `AbstractBoundaryCondition`s assigned to `Element`s
"
module StructuralBoundaryConditions

using Reexport
using Dictionaries: Dictionary, dictionary

using ..Utils
using ..Entities
using ..BoundaryConditions
using ..FixedFieldBoundaryConditions
using ..Nodes
using ..Meshes

const BCtoEntities{BC, E} = Dictionary{BC, Vector{E}} where {BC, E}

@reexport import ..BoundaryConditions: apply
@reexport import ..Entities: apply!
export StructuralBoundaryCondition, BCtoEntities, all_bcs, fixed_dof_bcs, load_bcs,
       displacement_bcs, element_bcs, face_bcs, node_bcs

"""
Structural boundary conditions.
A `StructuralBoundaryCondition` is a collection of `BoundaryConditions` defining the boundary conditions of the structure.
"""
Base.@kwdef struct StructuralBoundaryCondition{NB <: AbstractBoundaryCondition,
    FB <: AbstractBoundaryCondition,
    EB <: AbstractBoundaryCondition,
    N <: AbstractNode, F <: AbstractFace, E <: AbstractElement}
    "Maps each boundary conditions for a vector of nodes. "
    node_bcs::BCtoEntities{NB, N} = BCtoEntities{AbstractBoundaryCondition, AbstractNode}()
    "Maps each boundary conditions for a vector of faces."
    face_bcs::BCtoEntities{FB, F} = BCtoEntities{AbstractBoundaryCondition, AbstractFace}()
    "Maps each boundary conditions for a vector of elements. "
    element_bcs::BCtoEntities{EB, E} = BCtoEntities{
        AbstractBoundaryCondition, AbstractElement}()
end
function StructuralBoundaryCondition(pairs::Pair...)
    # Lenient constructor allowing to handle different concrete subtypes.
    StructuralBoundaryCondition([pairs...])
end
function StructuralBoundaryCondition(pairs::Vector{Pair{B,
        Vector{E}}}) where {B <: AbstractBoundaryCondition,
        E <: Union{AbstractEntity,
            AbstractNode}}
    node_bcs = BCtoEntities{AbstractBoundaryCondition, AbstractNode}()
    face_bcs = BCtoEntities{AbstractBoundaryCondition, AbstractFace}()
    element_bcs = BCtoEntities{AbstractBoundaryCondition, AbstractElement}()
    for (k, v) in pairs
        for vi in v
            if vi isa AbstractNode
                push!(get!(node_bcs, k, AbstractNode[]), vi)
            elseif vi isa AbstractFace
                push!(get!(face_bcs, k, AbstractFace[]), vi)
            elseif vi isa AbstractElement
                push!(get!(element_bcs, k, AbstractElement[]), vi)
            end
        end
    end
    StructuralBoundaryCondition(; node_bcs, face_bcs, element_bcs)
end

"Constructor for empty `StructuralBoundaryCondition` with a `Vector` of `AbstractBoundaryCondition`s `vbc`."
function StructuralBoundaryCondition(vbc::Vector{BC}) where {BC <:
                                                             AbstractBoundaryCondition}
    bcs_nodes = dictionary(map(bc -> bc => Vector{AbstractNode}(), vbc))
    bcs_faces = dictionary(map(bc -> bc => Vector{AbstractFace}(), vbc))
    bcs_elements = dictionary(map(bc -> bc => Vector{AbstractElement}(), vbc))
    StructuralBoundaryCondition(bcs_nodes, bcs_faces, bcs_elements)
end
function StructuralBoundaryCondition(bcs::AbstractBoundaryCondition...)
    StructuralBoundaryCondition(collect(bcs))
end

"Return the dictionary of node boundary conditions."
node_bcs(se::StructuralBoundaryCondition) = se.node_bcs

"Return the dictionary of face boundary conditions."
face_bcs(se::StructuralBoundaryCondition) = se.face_bcs

"Return the element boundary conditions."
element_bcs(se::StructuralBoundaryCondition) = se.element_bcs

"Return all the boundary conditions applied to the structure."
function all_bcs(se::StructuralBoundaryCondition)
    unique(vcat(collect(keys(node_bcs(se))), collect(keys(face_bcs(se))),
        collect(keys(element_bcs(se)))))
end

"Return the bc with the label `l` in the structural boundary conditions."
function Base.getindex(sb::StructuralBoundaryCondition, l::Label)
    first(filter(bc -> Symbol(label(bc)) == Symbol(l), all_bcs(sb)))
end

"Return a the entities where the boundary condition is applied."
function Base.getindex(sb::StructuralBoundaryCondition{NB, NF, EB},
        bc::BC) where
        {NB <: AbstractBoundaryCondition, NF <: AbstractBoundaryCondition,
        EB <: AbstractBoundaryCondition,
        BC <: AbstractBoundaryCondition}
    bc_entities = Vector{Union{AbstractNode, AbstractEntity}}()

    BC <: NB && bc ∈ keys(node_bcs(sb)) && push!(bc_entities, node_bcs(sb)[bc]...)
    BC <: NF && bc ∈ keys(face_bcs(sb)) && push!(bc_entities, face_bcs(sb)[bc]...)
    BC <: EB && bc ∈ keys(element_bcs(sb)) && push!(bc_entities, element_bcs(sb)[bc]...)

    isempty(bc_entities) ? throw(KeyError("Boundary condition $bc not found")) : bc_entities
end

function Base.push!(sb::StructuralBoundaryCondition, bc::AbstractBoundaryCondition,
        n::AbstractNode)
    push!(node_bcs(sb)[bc], n)
end
function Base.push!(sb::StructuralBoundaryCondition, bc::AbstractBoundaryCondition,
        e::AbstractElement)
    push!(element_bcs(sb)[bc], e)
end
function Base.push!(sb::StructuralBoundaryCondition, bc::AbstractBoundaryCondition,
        f::AbstractFace)
    push!(face_bcs(sb)[bc], f)
end

"Return a the boundary conditions with the label `l` in the structural boundary conditions."
function Base.getindex(sb::StructuralBoundaryCondition,
        l::L) where {L <: Union{Symbol, AbstractString}}
    bcs_label_l = collect(filter(bc -> Symbol(label(bc)) == Symbol(l), all_bcs(sb)))
    @assert length(bcs_label_l)==1 throw(ArgumentError("The label $l is not unique.
                                                          Please label each bc differently."))
    first(bcs_label_l)
end

"Return the boundary conditions applied to a node."
function Base.getindex(sb::StructuralBoundaryCondition, n::AbstractNode)
    keys(filter(x -> n ∈ x, node_bcs(sb)))
end

"Return the boundary conditions applied to a face."
function Base.getindex(sb::StructuralBoundaryCondition, f::AbstractFace)
    keys(filter(x -> f ∈ x, face_bcs(sb)))
end

"Return the boundary conditions applied to an element."
function Base.getindex(sb::StructuralBoundaryCondition, e::AbstractElement)
    keys(filter(x -> e ∈ x, element_bcs(sb)))
end

"Return a displacement boundary conditions."
function displacement_bcs(se::StructuralBoundaryCondition)
    vbc = Vector{AbstractDisplacementBoundaryCondition}()
    disp_bc = filter(bc -> bc isa AbstractDisplacementBoundaryCondition, all_bcs(se))
    push!(vbc, disp_bc...)
    unique(vbc)
end

"Return fixed dofs boundary conditions."
function fixed_dof_bcs(se::StructuralBoundaryCondition)
    vbc = Vector{FixedField}()
    fixed_bcs = filter(bc -> bc isa FixedField, all_bcs(se))
    push!(vbc, fixed_bcs...)
    unique!(vbc)
end

"Return the dofs that the boundary condition fixes."
function apply(bcs::StructuralBoundaryCondition, fbc::FixedField)
    # Extract nodes, faces and elements
    entities = bcs[fbc]
    dofs_to_delete = Dof[]
    for e in entities
        push!(dofs_to_delete, apply(fbc, e)...)
    end
    unique!(dofs_to_delete)
end
function apply(bcs::StructuralBoundaryCondition, f_bcs::Vector{<:FixedField})
    dofs_to_delete = Dof[]
    for fbc in f_bcs
        push!(dofs_to_delete, apply(bcs, fbc)...)
    end
    unique!(dofs_to_delete)
end

"Return the load boundary conditions."
function load_bcs(bcs::StructuralBoundaryCondition)
    load_bc = AbstractNeumannBoundaryCondition[]
    push!(load_bc, filter(bc -> bc isa AbstractNeumannBoundaryCondition, all_bcs(bcs))...)
    unique!(load_bc)
end

"Return a vector of dofs and load values applied by the boundary condition at time t."
function apply(
        bcs::StructuralBoundaryCondition, lbc::AbstractLoadBoundaryCondition, t::Real)
    # Extract nodes, faces and elements
    entities = bcs[lbc]
    dofs_to_load = Dof[]
    load_vec = Float64[]

    for e in entities
        dofs_lbc_e, load_vec_e = apply(lbc, e, t)
        push!(load_vec, load_vec_e...)
        push!(dofs_to_load, dofs_lbc_e...)
    end

    # Check if dofs are unique, if not add values
    num_dofs_to_load = length(dofs_to_load)
    unique_dofs_to_load = unique(dofs_to_load)
    if num_dofs_to_load == length(unique_dofs_to_load)
        dofs_to_load, load_vec
    else
        dofs_load_dict = dictionary(zip(unique_dofs_to_load, zeros(num_dofs_to_load)))
        for (dof, val) in zip(dofs_to_load, load_vec)
            dofs_load_dict[dof] += val
        end
        collect(keys(dofs_load_dict)), collect(values(dofs_load_dict))
    end
end

"Apply all the boundary conditions to the mesh."
function apply!(bcs::StructuralBoundaryCondition, m::AbstractMesh)
    apply_node_bcs!(bcs, m)
    apply_face_bcs!(bcs, m)
    apply_element_bcs!(bcs, m)
    bcs
end

function _assign_entities_to_bc!(bc_to_entities::BCtoEntities{BC},
        set_accessor::Function,
        mesh_entities::Vector{E},
        type_label::String,
        m::AbstractMesh) where {BC <: AbstractBoundaryCondition,
        E <: Union{AbstractEntity, AbstractNode}}
    sets = set_accessor(m)
    for (dbc, entities) in pairs(bc_to_entities)
        dbc_label = string(label(dbc))
        # Check if the boundary condition label is in the set
        # if not, the boundary condition is not applied
        if haskey(sets, dbc_label)
            for index in set_accessor(m, dbc_label)
                push!(entities, mesh_entities[index])
            end
        elseif !isempty(entities)
            @warn "The bc with label $dbc_label is already applied to the mesh."
        end
    end
end

"Apply node boundary conditions to the mesh. Assign entities with the same label in the
boundary conditions set."
function apply_node_bcs!(bcs::StructuralBoundaryCondition, m::AbstractMesh)
    _assign_entities_to_bc!(node_bcs(bcs),
        node_set,
        nodes(m),
        "node",
        m)

    _delete_empty_bc(node_bcs(bcs))
end
function _delete_empty_bc(bc_entity::BCtoEntities)
    for empty_bc in findall(isempty, bc_entity)
        delete!(bc_entity, empty_bc)
    end
end

"Apply face boundary conditions to the mesh. Assign entities with the same label in the
boundary conditions set."
function apply_face_bcs!(bcs::StructuralBoundaryCondition, m::AbstractMesh)
    _assign_entities_to_bc!(face_bcs(bcs),
        face_set,
        faces(m),
        "face",
        m)
    _delete_empty_bc(face_bcs(bcs))
end

"Apply element boundary conditions to the mesh. Assign entities with the same label in the
boundary conditions set."
function apply_element_bcs!(bcs::StructuralBoundaryCondition, m::AbstractMesh)
    _assign_entities_to_bc!(element_bcs(bcs),
        element_set,
        elements(m),
        "element",
        m)
    _delete_empty_bc(element_bcs(bcs))
end

"Replace the boundary condition with the same label as the new."
function Base.replace!(sb::StructuralBoundaryCondition,
        new_bc::AbstractBoundaryCondition,
        label::Label = label(new_bc))
    old_bc = sb[label]

    # node boundary conditions replacement
    node_boundary_conditions = node_bcs(sb)
    node_bc_type = eltype(keys(node_boundary_conditions))
    if (old_bc isa node_bc_type) && haskey(node_boundary_conditions, old_bc)
        old_nodes = node_boundary_conditions[old_bc]
        delete!(node_boundary_conditions, old_bc)
        insert!(node_boundary_conditions, new_bc, old_nodes)
    end

    # face boundary conditions replacement
    face_boundary_conditions = face_bcs(sb)
    face_bc_type = eltype(keys(face_boundary_conditions))
    if (old_bc isa face_bc_type) && haskey(face_boundary_conditions, old_bc)
        old_faces = face_boundary_conditions[old_bc]
        delete!(face_boundary_conditions, old_bc)
        insert!(face_boundary_conditions, new_bc, old_faces)
    end

    # element boundary conditions replacement
    element_boundary_conditions = element_bcs(sb)
    element_bc_type = eltype(keys(element_boundary_conditions))
    if (old_bc isa element_bc_type) && haskey(element_boundary_conditions, old_bc)
        old_elements = element_boundary_conditions[old_bc]
        delete!(element_boundary_conditions, old_bc)
        insert!(element_boundary_conditions, new_bc, old_elements)
    end
end

"Delete a boundary condition."
function Base.delete!(sb::StructuralBoundaryCondition,
        bc::AbstractBoundaryCondition)
    haskey(node_bcs(sb), bc) && delete!(node_bcs(sb), bc)
    haskey(element_bcs(sb), bc) && delete!(element_bcs(sb), bc)
    haskey(face_bcs(sb), bc) && delete!(face_bcs(sb), bc)
end
Base.delete!(sb::StructuralBoundaryCondition, label::Label) = delete!(sb, sb[label])

"Insert a boundary condition."
function Base.insert!(sb::StructuralBoundaryCondition,
        bc::AbstractBoundaryCondition,
        bc_entities::Vector{<:AbstractNode})
    insert!(node_bcs(sb), bc, bc_entities)
end
function Base.insert!(sb::StructuralBoundaryCondition,
        bc::AbstractBoundaryCondition,
        bc_entities::AbstractNode...)
    insert!(sb, bc, collect(bc_entities))
end
function Base.insert!(sb::StructuralBoundaryCondition,
        bc::AbstractBoundaryCondition,
        bc_entities::Vector{<:AbstractFace})
    insert!(face_bcs(sb), bc, bc_entities)
end
function Base.insert!(sb::StructuralBoundaryCondition,
        bc::AbstractBoundaryCondition,
        bc_entities::AbstractFace...)
    insert!(sb, bc, collect(bc_entities))
end
function Base.insert!(sb::StructuralBoundaryCondition,
        bc::AbstractBoundaryCondition,
        bc_entities::Vector{<:AbstractElement})
    insert!(element_bcs(sb), bc, bc_entities)
end
function Base.insert!(sb::StructuralBoundaryCondition,
        bc::AbstractBoundaryCondition,
        bc_entities::AbstractElement...)
    insert!(sb, bc, collect(bc_entities))
end

end # module
