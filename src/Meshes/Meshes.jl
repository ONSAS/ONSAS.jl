"""
Module defining a mesh interface, `AbstractMesh`, and the default implementation `Mesh`.

Each mesh consists of different kinds of entities, stores as dense arrays:

- Nodes.
- Entities.
- Faces.

In addition to the entities, metadata mapping names to entity indices (`EntitySet`) can be provided.
"""
module Meshes

using Dictionaries
using Reexport

using ..Entities
using ..Nodes
using ..Utils

@reexport import ..Entities: dimension, dofs, nodes, num_nodes
@reexport import ..Nodes: set_dofs!

export AbstractMesh, Mesh, EntitySet, connectivity, faces, face_set, node, element,
       elements,
       element_set, num_dofs, num_elements, node_set, add_node_to_set!, add_element_to_set!,
       add_face_to_set!, add_entity_to_set!, node_matrix

"""
Abstract mesh of dimension `dim`.

**Common methods**
* [`connectivity`](@ref)
* [`dimension`](@ref)
* [`dofs`](@ref)
* [`num_dofs`](@ref)
* [`set_dofs`](@ref)
* [`element`](@ref)
* [`elements`](@ref)
* [`face`](@ref)
* [`faces`](@ref)
* [`num_elements`](@ref)
* [`node`](@ref)
* [`nodes`](@ref)
* [`num_nodes`](@ref)
* [`node_matrix`](@ref)
"""
abstract type AbstractMesh{dim} end

"Return the dimension of an `AbstractMesh`."
dimension(::AbstractMesh{dim}) where {dim} = dim

"""
Return an iterator over the degrees of freedom for each node.
The `i`-th entry contains the `Dof`s of node with index `i` in the mesh.
"""
dofs(m::AbstractMesh) = Iterators.map(dofs, nodes(m))

"""
Return the number of dofs defined in the mesh.
This function assumes that `Dof`s indexes start from `Dof(1)`.
"""
function num_dofs(m::AbstractMesh)::Int
    max_dof = 0
    for n in nodes(m)
        # For each node, obtain the maximum dof over each field, then reduce over all fields.
        max_n = mapreduce(maximum, max, dofs(n); init = 0)
        max_dof = max(max_dof, max_n)
    end
    max_dof
end
function num_dofs(m::AbstractMesh, f::Field)::Int
    dof_indices = Int[]

    for n in nodes(m)
        # Append the DOF indices for the specified field `f` for each node
        append!(dof_indices, dofs(n, f))
    end

    # Use `unique` to remove duplicate DOF indices, then count them
    return length(unique!(dof_indices))
end

"Set `dofs_per_node` degrees of freedom per node with the given symbol to all nodes of the mesh."
function set_dofs!(m::AbstractMesh, dof_symbol::Field, dofs_per_node::Int)
    # Check that there are no nodes with a dof associated to the given symbol.
    dof_exists = any(n -> dof_symbol ∈ keys(dofs(n)), nodes(m))
    if dof_exists
        throw(ArgumentError("Dof symbol $dof_symbol already exists."))
    end

    max_dof = num_dofs(m)
    for (i, n) in enumerate(nodes(m))
        h = max_dof + i * dofs_per_node
        node_dofs = [Dof(k) for k in (1 + h - dofs_per_node):h]
        set_dofs!(n, dof_symbol, node_dofs)
    end
end

"Return the element at index `i`."
node(m::AbstractMesh, i::Int) = nodes(m)[i]

"Return the array of nodes."
nodes(m::AbstractMesh) = m.nodes

"Return the number of nodes."
num_nodes(m::AbstractMesh) = length(nodes(m))

"Return the array of faces."
faces(m::AbstractMesh) = m.faces

"Return the number of faces."
num_faces(m::AbstractMesh) = length(faces(m))

"Return the array of elements."
elements(m::AbstractMesh) = m.elements

"Return the element at index `i`."
element(m::AbstractMesh, i::Int) = elements(m)[i]

"Return the number of elements."
num_elements(m::AbstractMesh) = length(elements(m))

"""
    push!(m::AbstractMesh, n::AbstractNode)
    push!(m::AbstractMesh, e::AbstractElement)
    push!(m::AbstractMesh, f::AbstractFace)

Push a node, element or face to the mesh.
Important: these methods do *not* update the nodes array in the mesh.
Used mostly for programmatic construction of meshes via GMSH.
"""
Base.push!(m::AbstractMesh, n::AbstractNode) = push!(nodes(m), n)
Base.push!(m::AbstractMesh, e::AbstractElement) = push!(elements(m), e)
Base.push!(m::AbstractMesh, f::AbstractFace) = push!(faces(m), f)

# TODO Dispatch on type.
Base.length(m::AbstractMesh, ::AbstractNode) = length(nodes(m))
Base.length(m::AbstractMesh, ::AbstractElement) = length(elements(m))
Base.length(m::AbstractMesh, ::AbstractFace) = length(faces(m))

"Used to designate node, element and face sets mapping an entity string to a set of indices."
const EntitySet = Dictionary{String, Set{Int}}

"""
A `Mesh` is a collection of `Element`s, `Face`s and `Node`s that cover the discretized domain,
together with Sets of elements and nodes.

### Methods

* [`node_set`](@ref)
* [`add_node_to_set!`](@ref)
* [`element_set`](@ref)
* [`add_element_to_set!`](@ref)
* [`face_set`](@ref)
* [`add_face_to_set!`](@ref)
"""
struct Mesh{dim, N <: AbstractNode{dim}, E <: AbstractElement, F <: AbstractFace, EX} <:
       AbstractMesh{dim}
    "`Node`s of the mesh with dimension `dim`."
    nodes::Vector{N}
    "`Element`s of the mesh."
    elements::Vector{E}
    "`Face`s of the mesh."
    faces::Vector{F}
    "Sets of `Node`s mapped to a string key."
    node_sets::EntitySet
    "Sets of `Element`s mapped to a string key."
    element_sets::EntitySet
    "Sets of `Face`s mapped to a string key."
    face_sets::EntitySet
    "Additional data or mesh info."
    extra::EX
end
function Mesh(; nodes::Vector{N} = Vector{AbstractNode}(),
        elements::Vector{E} = Vector{AbstractElement}(),
        faces::Vector{F} = Vector{AbstractFace}(),
        node_sets::EntitySet = EntitySet(),
        face_sets::EntitySet = EntitySet(),
        element_sets::EntitySet = EntitySet(),
        extra::EX = nothing) where
        {N <: AbstractNode, F <: AbstractFace, E <: AbstractElement, EX}
    Mesh(nodes, elements, faces, node_sets, element_sets, face_sets, extra)
end

"Return the mesh node coordinates matrix. Each row is a node, each column a coordinate."
function node_matrix(mesh::Mesh{dim, T}) where {dim, T}
    nodes_coords_matrix = Matrix{eltype(T)}(undef, (dim, num_nodes(mesh)))
    for (i, n) in enumerate(nodes(mesh))
        nodes_coords_matrix[:, i] = coordinates(n)
    end
    nodes_coords_matrix
end

"Return the mesh connectivity."
function connectivity(mesh::Mesh{dim, T}) where {dim, T}

    # Check if a already contains the connectivity
    hasproperty(mesh.extra, :connectivity) && return mesh.extra.connectivity

    connectivity = Vector{Vector{Int}}(undef, num_elements(mesh))

    enumerate_nodes = Dictionary{AbstractNode, Int}()
    for (i, n) in enumerate(nodes(mesh))
        get!(enumerate_nodes, n, i)
    end

    for (i, e) in enumerate(elements(mesh))
        connectivity[i] = [enumerate_nodes[n] for n in nodes(e)]
    end
    connectivity
end

"Return the set of named nodes."
node_set(m::Mesh) = m.node_sets

"Return the set of nodes associated to a given name."
node_set(m::Mesh, name::String) = get(node_set(m), name, nothing)

"Return a `Vector` of `Node`s  with `node_set_name` in the `Mesh` `m`."
function nodes(m::Mesh, name::String)
    idx = node_set(m, name)
    mesh_nodes = nodes(m)
    # TODO vcat ?
    [mesh_nodes[i] for i in idx]
end

"Register a node's name in the mesh."
function add_node_to_set!(m::Mesh, name::String, n::AbstractNode)
    idx = findfirst(==(n), nodes(m))
    isnothing(idx) && error("Node $n not found in the mesh.")
    add_node_to_set!(node_set(m), name, idx)
end
function add_node_to_set!(m::Mesh, name::String, idx::Int)
    add_node_to_set!(node_set(m), name, idx)
end
function add_node_to_set!(node_sets::EntitySet, name::AbstractString, idx::Int)
    if haskey(node_sets, name)
        push!(node_sets[name], idx)
    else
        insert!(node_sets, name, Set(idx))
    end
    node_sets
end

"Return the set of named elements."
element_set(m::Mesh) = m.element_sets

"Return the set of elements associated to a given name."
element_set(m::Mesh, name::String) = get(element_set(m), name, nothing)

"Return a `Vector` of `Element`s  with `element_set_name` in the `Mesh` `m`."
function elements(m::Mesh, name::String)
    idx = element_set(m, name)
    elems = elements(m)
    [elems[i] for i in idx]
end

function Base.show(io::IO, m::Mesh)
    nnodes = num_nodes(m)
    nelems = num_elements(m)
    nfaces = num_faces(m)
    println("• Mesh with $nnodes nodes, $nelems elements and $nfaces faces.")
end

"Register an element's name in the mesh."
function add_element_to_set!(m::Mesh, name::AbstractString, e::AbstractElement)
    idx = findfirst(==(e), elements(m))
    isnothing(idx) && error("Element $e not found in the mesh.")
    add_element_to_set!(element_set(m), name, idx)
end
function add_element_to_set!(m::Mesh, name::AbstractString, idx::Int)
    add_element_to_set!(element_set(m), name, idx)
end
function add_element_to_set!(element_sets::EntitySet, name::AbstractString, idx::Int)
    if haskey(element_sets, name)
        push!(element_sets[name], idx)
    else
        insert!(element_sets, name, Set(idx))
    end
    element_sets
end

"Return the set of named faces."
face_set(m::Mesh) = m.face_sets

"Return the set of faces associated to a given name."
face_set(m::Mesh, name::String) = get(face_set(m), name, nothing)

"Return a `Vector` of `Face`s  with `name` in the `Mesh` `m`."
function faces(m::Mesh, name::String)
    face_indexes = face_set(m, name)
    mesh_faces = faces(m)
    [mesh_faces[i] for i in face_indexes]
end

"Register a face's name in the mesh."
function add_face_to_set!(m::Mesh, name::AbstractString, f::AbstractFace)
    idx = findfirst(==(f), faces(m))
    isnothing(idx) && error("Face $f not found in mesh.")
    add_face_to_set!(face_set(m), name, idx)
end
function add_face_to_set!(m::Mesh, name::AbstractString, idx::Int)
    add_face_to_set!(face_set(m), name, idx)
end
function add_face_to_set!(face_sets::EntitySet, name::AbstractString, idx::Int)
    if haskey(face_sets, name)
        push!(face_sets[name], idx)
    else
        insert!(face_sets, name, Set(idx))
    end
    face_sets
end

"Add an entity to a set, dispatching on the entity type."
function add_entity_to_set!(mesh::AbstractMesh, entity_type_label::AbstractString,
        entity_position::Int,
        ::AbstractNode)
    add_node_to_set!(mesh, entity_type_label, entity_position)
end
function add_entity_to_set!(mesh::AbstractMesh, entity_type_label::AbstractString,
        entity_position::Int,
        ::AbstractFace)
    add_face_to_set!(mesh, entity_type_label, entity_position)
end
function add_entity_to_set!(mesh::AbstractMesh, entity_type_label::AbstractString,
        entity_position::Int,
        ::AbstractElement)
    add_element_to_set!(mesh, entity_type_label, entity_position)
end

"Replace a node in the mesh."
function Base.replace!(mesh::AbstractMesh{dim},
        node_idx::Int,
        node_coordinates::Point{dim}) where {dim}
    mesh_nodes = nodes(mesh)
    old_node = mesh_nodes[node_idx]
    new_node = Node(node_coordinates, dofs(old_node))
    mesh_nodes[node_idx] = new_node
    # Check faces and elements nodes are views from the mesh nodes.
    if !isempty(faces(mesh))
        @assert nodes(rand(faces(mesh))) isa SubArray
    end

    if !isempty(elements(mesh))
        @assert nodes(rand(elements(mesh))) isa SubArray
    end
    new_node
end

end # module
