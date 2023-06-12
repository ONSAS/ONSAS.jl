module Structures

using Reexport
using Dictionaries: Dictionary, dictionary

using ..Entities
using ..Meshes
using ..Gmsh
using ..BoundaryConditions
using ..StructuralModel
using ..StructuralMaterials
using ..StructuralEntities
using ..StructuralBoundaryConditions
using ..Utils
using ..Nodes

@reexport import ..Handlers: PointEvalHandler

export Structure

"""
An `Structure` object facilitates the process of assembling and creating the structural analysis. 
"""
struct Structure{dim,MESH,MAT,E,NB,LB} <: AbstractStructure{dim,MAT,E}
    "Stores the structure's mesh."
    mesh::MESH
    "Stores the structure's materials and elements assignments."
    materials::StructuralMaterial{MAT,E}
    "Stores the structure's boundary conditions and elements assignments."
    bcs::StructuralBoundaryCondition{NB,LB}
    "Stores the structure's free degrees of freedom."
    free_dofs::Vector{Dof}
    function Structure(mesh::MESH,
                       materials::StructuralMaterial{MAT,E},
                       bcs::StructuralBoundaryCondition{NB,LB},
                       free_dofs::Vector{Dof}) where {dim,MESH<:AbstractMesh{dim},MAT,E,NB,LB}
        new{dim,MESH,MAT,E,NB,LB}(mesh, materials, bcs, free_dofs)
    end
end

"Constructor with  `StructuralMaterial` `materials`,  `StructuralBoundaryCondition` `bcs` 
and `AbstractMesh` `mesh` seting fixed dofs with `FixedDofBoundaryCondition` defined in `bcs`"
function Structure(mesh::AbstractMesh{dim},
                   materials::StructuralMaterial{M,E},
                   bcs::StructuralBoundaryCondition{NB,LB}) where {dim,M,E,NB,LB}
    default_free_dofs = Vector{Dof}()
    for node_dofs in dofs(mesh)
        [push!(default_free_dofs, vec_dof...) for vec_dof in collect(values(node_dofs))]
    end

    fixed_dofs = apply(bcs, fixed_dof_bcs(bcs))

    deleteat!(default_free_dofs, findall(x -> x in fixed_dofs, default_free_dofs))

    Structure(mesh, materials, bcs, default_free_dofs)
end

"Constructor of a `Structure` given a `MshFile` `msh_file`, `StructuralMaterial` `materials`,  `StructuralBoundaryCondition` `bcs`."
function Structure(msh_file::MshFile,
                   materials::StructuralMaterial, bcs::StructuralBoundaryCondition,
                   s_entities::StructuralEntity,
                   dofs_to_dim::Dictionary{Field,<:Integer}=dictionary([:u => 3]))
    nodes = msh_file.vec_nodes
    mesh = Mesh(; nodes)

    # Loop over all physical entities
    for (entity_index, entity_nodes_indexes) in enumerate(msh_file.connectivity)

        # Create entity and push it into the mesh
        nodes_entity = view(nodes, entity_nodes_indexes)
        entity_type_label = entity_label(msh_file, entity_index)
        # Check if the entity is a node, if not add it to the mesh
        if entity_type_label == "node"
            entity = nodes_entity[]
        else
            entity_type = s_entities[entity_type_label]
            entity = create_entity(entity_type, nodes_entity)
            push!(mesh, entity)
        end

        # Find material and push 
        material_type_label = material_label(msh_file, entity_index)
        # If has material defined is an element not a surface
        if ~isempty(material_type_label)
            material_type = materials[material_type_label]
            push!(materials[material_type], entity)
        end

        # Find boundary conditions 
        bc_type_label = bc_label(msh_file, entity_index)
        if !isempty(bc_type_label)
            bc_type = bcs[bc_type_label]
            push!(bcs, bc_type, entity)
        end
    end

    for (dof_symbol, dof_dim) in pairs(dofs_to_dim)
        apply!(mesh, dof_symbol, dof_dim)
    end

    Structure(mesh, materials, bcs)
end

function Base.show(io::IO, s::Structure)
    ndofs = length(s.free_dofs)
    println("• Structure with $ndofs free dofs.")
    show(io, s.mesh)
end

"Constructor of a `PointEvalHandler` from a `Structure` and a `AbstractVector` of `Point`s ."
function PointEvalHandler(s::Structure, vec_points::AbstractVector{P}) where {T,P<:Point{T}}
    PointEvalHandler(mesh(s), vec_points)
end

"Replace the `AbstractMaterial` with the label `l` for the new material `new_m` in the `Structure`.
The previous material `Element`s are assigned to the new."
function Base.replace!(s::Structure,
                       new_material::AbstractMaterial,
                       label::Label=label(new_material))
    replace!(materials(s), new_material, label)
end

end # module