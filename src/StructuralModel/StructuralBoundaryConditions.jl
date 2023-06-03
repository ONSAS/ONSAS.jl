using Reexport

using ..Utils
using ..Entities
using ..BoundaryConditions
using ..FixedDofBoundaryConditions
using ..Nodes

@reexport import ..BoundaryConditions: apply
@reexport import ..Entities: apply!

export StructuralBoundaryCondition, all_bcs, node_bcs, face_bcs, element_bcs, displacement_bcs,
       load_bcs, fixed_dof_bcs

""" Structural boundary conditions.
A `StructuralBoundaryCondition` is a collection of `BoundaryConditions` defining the boundary conditions of the structure.
"""
Base.@kwdef struct StructuralBoundaryCondition{NB<:AbstractBoundaryCondition,
                                               FB<:AbstractBoundaryCondition,
                                               EB<:AbstractBoundaryCondition,
                                               N<:AbstractNode,F<:AbstractFace,E<:AbstractElement}
    "Maps each boundary conditions for a vector of nodes. "
    node_bcs::Dictionary{NB,Vector{N}} = Dictionary{AbstractBoundaryCondition,Vector{AbstractNode}}()
    "Maps each boundary conditions for a vector of faces."
    face_bcs::Dictionary{FB,Vector{F}} = Dictionary{AbstractBoundaryCondition,Vector{AbstractFace}}()
    "Maps each boundary conditions for a vector of elements. "
    element_bcs::Dictionary{EB,Vector{E}} = Dictionary{AbstractBoundaryCondition,
                                                       Vector{AbstractElement}}()
end

"Constructor for empty `StructuralBoundaryCondition` with a `Vector` of `AbstractBoundaryCondition`s `vbc`."
function StructuralBoundaryCondition(vbc::Vector{BC}) where {BC<:AbstractBoundaryCondition}
    bcs_nodes = dictionary(map(bc -> bc => Vector{AbstractNode}(), vbc))
    bcs_faces = dictionary(map(bc -> bc => Vector{AbstractFace}(), vbc))
    bcs_elements = dictionary(map(bc -> bc => Vector{AbstractElement}(), vbc))
    StructuralBoundaryCondition(bcs_nodes, bcs_faces, bcs_elements)
end
function StructuralBoundaryCondition(bcs::AbstractBoundaryCondition...)
    StructuralBoundaryCondition(collect(bcs))
end

"Return the `BoundaryCondition` with the label `l` in the `StructuralBoundaryCondition` `sb`."
function Base.getindex(sb::StructuralBoundaryCondition, l::Label)
    first(filter(bc -> Symbol(label(bc)) == Symbol(l), all_bcs(sb)))
end

"Return the `Vector` of `Node`s and `Element`s where the `BoundaryCondition` `bc` is imposed."
function Base.getindex(sb::StructuralBoundaryCondition{NB,NF,EB},
                       bc::BC) where
         {NB<:AbstractBoundaryCondition,NF<:AbstractBoundaryCondition,EB<:AbstractBoundaryCondition,
          BC<:AbstractBoundaryCondition}
    bc_entities = Vector{Union{AbstractNode,AbstractEntity}}()

    BC <: NB && bc ∈ keys(node_bcs(sb)) && push!(bc_entities, node_bcs(sb)[bc]...)
    BC <: NF && bc ∈ keys(face_bcs(sb)) && push!(bc_entities, face_bcs(sb)[bc]...)
    BC <: EB && bc ∈ keys(element_bcs(sb)) && push!(bc_entities, element_bcs(sb)[bc]...)

    isempty(bc_entities) ? throw(KeyError("Boundary condition $bc not found")) : bc_entities
end

"Return a the `BoundaryConditions` with the label `l` in the `StructuralEntity` `sb`."
function Base.getindex(sb::StructuralBoundaryCondition,
                       l::L) where {L<:Union{Symbol,AbstractString}}
    bcs_label_l = collect(filter(bc -> Symbol(label(bc)) == Symbol(l), all_bcs(sb)))
    @assert length(bcs_label_l) == 1 throw(ArgumentError("The label $l is not unique. Please label each bc differently."))
    first(bcs_label_l)
end

"Return the `Vector` of `BoundaryConditions`s applied to the `AbstractNode` `n`."
function Base.getindex(sb::StructuralBoundaryCondition, n::AbstractNode)
    keys(filter(x -> n ∈ x, node_bcs(sb)))
end

"Return the `Vector` of `BoundaryConditions`s applied to the `AbstractFace` `f`."
function Base.getindex(sb::StructuralBoundaryCondition, f::AbstractFace)
    keys(filter(x -> f ∈ x, face_bcs(sb)))
end

"Return the `Vector` of `BoundaryConditions`s applied to the `AbstractElement` `e`."
function Base.getindex(sb::StructuralBoundaryCondition, e::AbstractElement)
    keys(filter(x -> e ∈ x, element_bcs(sb)))
end

"Pushes to the `AbstractNode` `n` the `BoundaryCondition` `bc` in the `StructuralBoundaryCondition` `sb`."
function Base.push!(sb::StructuralBoundaryCondition, bc::AbstractBoundaryCondition,
                    n::AbstractNode)
    push!(node_bcs(sb)[bc], n)
end

"Pushes to the `AbstractElement` `e` the `BoundaryCondition` `bc` in the `StructuralBoundaryCondition` `sb`."
function Base.push!(sb::StructuralBoundaryCondition, bc::AbstractBoundaryCondition,
                    e::AbstractElement)
    push!(element_bcs(sb)[bc], e)
end

"Pushes to the `AbstractFace` `f` the `BoundaryCondition` `bc` in the `StructuralBoundaryCondition` `sb`."
function Base.push!(sb::StructuralBoundaryCondition, bc::AbstractBoundaryCondition,
                    f::AbstractFace)
    push!(face_bcs(sb)[bc], f)
end

"Return the dictionary of `BoundaryConditions`s applied to `Node`s."
node_bcs(se::StructuralBoundaryCondition) = se.node_bcs

"Return the dictionary of `BoundaryConditions`s applied to `Face`s."
face_bcs(se::StructuralBoundaryCondition) = se.face_bcs

"Return the `Dictionary` of `BoundaryConditions`s applied to `Element`s."
element_bcs(se::StructuralBoundaryCondition) = se.element_bcs

"Return all `BoundaryConditions`s defined into `StructuralBoundaryCondition`."
function all_bcs(se::StructuralBoundaryCondition)
    unique(vcat(collect(keys(node_bcs(se))), collect(keys(face_bcs(se))),
                collect(keys(element_bcs(se)))))
end

"Return a `Vector` of `AbstractDisplacementBoundaryCondition`s applied to `Node`s and `Element`s in the `StructuralBoundaryCondition` `se`."
function displacement_bcs(se::StructuralBoundaryCondition)
    vbc = Vector{AbstractDisplacementBoundaryCondition}()
    disp_bc = filter(bc -> bc isa AbstractDisplacementBoundaryCondition, all_bcs(se))
    push!(vbc, disp_bc...)
    unique(vbc)
end

"Return a `Vector` of `FixedDof`s applied to `Node`s, `Face`s and `Element`s. in the `StructuralBoundaryCondition` `se`."
function fixed_dof_bcs(se::StructuralBoundaryCondition)
    vbc = Vector{FixedDof}()
    fixed_bcs = filter(bc -> bc isa FixedDof, all_bcs(se))
    push!(vbc, fixed_bcs...)
    unique(vbc)
end

"Return a `Vector` of `Dof`s to delete given a `FixedDof` and a set of `StructuralBoundaryCondition` `bcs`."
function apply(bcs::StructuralBoundaryCondition, fbc::FixedDof)
    # Extract nodes, faces and elements 
    entities = bcs[fbc]
    dofs_to_delete = Dof[]
    [push!(dofs_to_delete, apply(fbc, e)...) for e in entities]
    unique!(dofs_to_delete)
end

"Return a `Vector` of `FixedDof` `f_bcs` and a set of `StructuralBoundaryCondition` `bcs`."
function apply(bcs::StructuralBoundaryCondition, f_bcs::Vector{<:FixedDof})
    dofs_to_delete = Dof[]
    [push!(dofs_to_delete, apply(bcs, fbc)...) for fbc in f_bcs]
    unique!(dofs_to_delete)
end

"Return a `Vector` of `FixedDof`s applied to `Node`s and `Element`s in the `StructuralBoundaryCondition` `bcs`."
function load_bcs(bcs::StructuralBoundaryCondition)
    vbc = Vector{AbstractLoadBoundaryCondition}()
    load_bcs = filter(bc -> bc isa AbstractLoadBoundaryCondition, all_bcs(bcs))
    push!(vbc, load_bcs...)
    unique!(vbc)
end

"Return a `Vector` of `Dof`s and values to apply a `LoadDofBoundaryCondition` and a set of `StructuralBoundaryCondition` `bcs`
at time `t`."
function apply(bcs::StructuralBoundaryCondition, lbc::AbstractLoadBoundaryCondition, t::Real)
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
        [dofs_load_dict[dof] += val for (dof, val) in zip(dofs_to_load, load_vec)]
        collect(keys(dofs_load_dict)), collect(values(dofs_load_dict))
    end
end

"Apply the `StructuralBoundaryCondition` to the `AbstractMesh` `m`. For this is required sets 
into the `Mesh` and the corresponding boundary condition labels declared in `bcs`."
function apply!(bcs::StructuralBoundaryCondition, m::AbstractMesh)
    apply_node_bcs!(bcs, m)
    apply_face_bcs!(bcs, m)
    apply_element_bcs!(bcs, m)
    _delete_empty_bcs!(bcs)
    bcs
end

"Apply `Node` boundary conditions given the `AbstractMesh` `m` to the `StructuralBoundaryCondition` `bcs`."
function apply_node_bcs!(bcs::StructuralBoundaryCondition, m::AbstractMesh)
    # Assign entities to the node boundary conditions
    vec_nodes = nodes(m)
    node_sets = node_set(m)
    node_boundary_conditions = node_bcs(bcs)

    for (dbc, entities) in pairs(node_boundary_conditions)
        dbc_label = string(label(dbc))
        # Check if the boundary conditions label is the mesh node set 
        # if not, the boundary condition is not applied
        if haskey(node_sets, dbc_label)
            [push!(entities, vec_nodes[node_index]) for node_index in node_set(m, dbc_label)]
        end
    end
end

"Apply `Face` boundary conditions given the `AbstractMesh` `m` to the `StructuralBoundaryCondition` `bcs`."
function apply_face_bcs!(bcs::StructuralBoundaryCondition, m::AbstractMesh)
    # Assign entities to the node boundary conditions
    vec_faces = faces(m)
    face_sets = face_set(m)
    face_boundary_conditions = face_bcs(bcs)

    for (dbc, entities) in pairs(face_boundary_conditions)
        dbc_label = string(label(dbc))
        # Check if the boundary conditions label is the mesh face set 
        # if not, the boundary condition is not applied
        if haskey(face_sets, dbc_label)
            [push!(entities, vec_faces[face_index]) for face_index in face_set(m, dbc_label)]
        end
    end
end

"Apply `Element` boundary conditions given the `AbstractMesh` `m` to the `StructuralBoundaryCondition` `bcs`."
function apply_element_bcs!(bcs::StructuralBoundaryCondition, m::AbstractMesh)
    # Assign entities to the node boundary conditions
    vec_elements = elements(m)
    element_sets = element_set(m)
    element_boundary_conditions = element_bcs(bcs)

    for (dbc, entities) in pairs(element_boundary_conditions)
        dbc_label = string(label(dbc))
        # Check if the boundary conditions label is the mesh element set 
        # if not, the boundary condition is not applied
        if haskey(element_sets, dbc_label)
            [push!(entities, vec_elements[element_index])
             for element_index in element_set(m, dbc_label)]
        end
    end
end

"Delete empty boundary conditions from the `StructuralBoundaryCondition` `bcs`."
function _delete_empty_bcs!(bcs::StructuralBoundaryCondition)
    bcs_dicts = [node_bcs(bcs), face_bcs(bcs), element_bcs(bcs)]
    for bc_dict in bcs_dicts
        for (dbc, entities) in pairs(bc_dict)
            isempty(entities) && delete!(bc_dict, dbc)
        end
    end
    bcs
end
