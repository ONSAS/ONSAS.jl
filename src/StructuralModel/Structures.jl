"""
Module defining structure entities interface.
Each Structure consists of a data type with materials, boundary conditions, mesh and a vector of free dofs.
"""
module Structures

using Reexport
using Dictionaries: Dictionary, dictionary

using ..Entities
using ..Nodes
using ..Meshes
using ..Gmsh
using ..BoundaryConditions
using ..Materials
using ..StructuralMaterials
using ..StructuralEntities
using ..StructuralBoundaryConditions
using ..Utils

@reexport import ..Meshes: dofs, num_dofs, nodes, num_nodes, faces, num_faces,
                           elements, num_elements
@reexport import ..Handlers: PointEvalHandler, mesh

export AbstractStructure, Structure, materials, boundary_conditions, num_free_dofs, free_dofs

# ==========
# Structure
# ==========

"""
Abstract supertype to define a new structure.
An `AbstractStructure` object facilitates the process of assigning materials, elements,
and boundary conditions to the mesh. Moreover this struct is used to solve an specific
structural analysis.

**Abstract Methods**

## Mesh:
* [`dofs`](@ref)
* [`num_dofs`](@ref)
* [`free_dofs`](@ref)
* [`num_free_dofs`](@ref)
* [`nodes`](@ref)
* [`num_nodes`](@ref)
* [`faces`](@ref)
* [`num_faces`](@ref)
* [`elements`](@ref)
* [`num_elements`](@ref)
* [`mesh`](@ref)

## Boundary conditions:
* [`boundary_conditions`](@ref)
* [`displacement_bcs`](@ref)
* [`load_bcs`](@ref)
* [`node_bcs`](@ref)
* [`element_bcs`](@ref)

## Materials:
* [`materials`](@ref)

"""
abstract type AbstractStructure{dim,M,E} end

# Mesh methods
"Return the structure mesh."
mesh(s::AbstractStructure) = s.mesh

"Return the structure dofs."
dofs(s::AbstractStructure) = dofs(mesh(s))

"Return structure dofs."
num_dofs(s::AbstractStructure) = num_dofs(mesh(s))

"Return the structure free dofs."
free_dofs(s::AbstractStructure) = s.free_dofs

"Return the number of free dofs of the structure."
num_free_dofs(s::AbstractStructure) = length(free_dofs(s))

"Return the structure nodes."
nodes(s::AbstractStructure) = nodes(mesh(s))

"Return the structure number of nodes."
num_nodes(s::AbstractStructure) = num_nodes(mesh(s))

"Return the structure elements."
elements(s::AbstractStructure) = elements(mesh(s))

"Return the structure number of elements."
num_elements(s::AbstractStructure) = num_elements(mesh(s))

"Return the structure faces."
faces(s::AbstractStructure) = faces(mesh(s))

"Return the number of `Face`s of the `AbstractStructure` `s`"
num_faces(s::AbstractStructure) = num_faces(mesh(s))

# Boundary Conditions
"Return the structure boundary conditions."
boundary_conditions(s::AbstractStructure) = s.bcs

# Materials
"Return the structure materials."
materials(s::AbstractStructure) = s.materials

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

"Constructor for a structure with mesh, materials and boundary conditions."
function Structure(mesh::AbstractMesh{dim},
                   materials::StructuralMaterial{M,E},
                   bcs::StructuralBoundaryCondition{NB,LB}) where {dim,M,E,NB,LB}
    # Retrieve all (free) dofs in the mesh.
    default_free_dofs = Vector{Dof}()
    for n in nodes(mesh)
        append!(default_free_dofs, reduce(vcat, dofs(n)))
    end
    # Obtain the dofs that shall be fixed according to the given boundary condition.
    fixed_dofs = apply(bcs, fixed_dof_bcs(bcs))
    # Remove the fixed dofs form the free dofs array.
    setdiff!(default_free_dofs, fixed_dofs)
    Structure(mesh, materials, bcs, default_free_dofs)
end

"Constructor for a structure given a mesh file, and their corresponding structural materials,
structural elements and structural boundary conditions."
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
            material_type = materials[Symbol(material_type_label)]
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
        set_dofs!(mesh, dof_symbol, dof_dim)
    end

    Structure(mesh, materials, bcs)
end

function Base.show(io::IO, s::Structure)
    ndofs = num_free_dofs(s)
    println("• Structure with $ndofs free dofs.")
    show(io, s.mesh)
end

"Constructor of a point eval handler from a structure and a vector of points ."
function PointEvalHandler(s::Structure, vec_points::AbstractVector{P}) where {T,P<:Point{T}}
    PointEvalHandler(mesh(s), vec_points)
end

"Replace a material with the a label for the new material in the structure.
The previous material element's are assigned to the new."
function Base.replace!(s::Structure,
                       new_material::AbstractMaterial,
                       label::Label=label(new_material))
    replace!(materials(s), new_material, label)
end

end # module
