"""
Module defining meshes entities interface.
Each mesh consists of a data type with the nodes and elements. Moreover, different sets of nodes and elements can be defined.
"""
module Meshes

using Dictionaries
using Reexport

@reexport using ..Elements
@reexport import ..Elements: apply!, dimension, dofs, nodes, num_nodes

export AbstractMesh, Mesh, faces, face_set, element, elements, element_set, num_dofs, num_elements,
       node_set, add_node_to_set!, add_element_to_set!, add_face_to_set!

""" Abstract supertype for all meshes.

The following methods are provided by the interface:


**Common methods:**

* [`dimension`](@ref)
* [`dofs`](@ref)
* [`num_dofs`](@ref)
* [`elements`](@ref)
* [`faces`](@ref)
* [`num_elements`](@ref)
* [`nodes`](@ref)
* [`num_nodes`](@ref)
"""
abstract type AbstractMesh{dim} end

"Return the dimension of an `AbstractMesh`."
dimension(::AbstractMesh{dim}) where {dim} = dim

"Return the `AbstractMesh` vector of `Dof`s. 
Entry `i` contains the `Dof`s of node with index `i` in the `AbstractMesh` vector of nodes."
dofs(m::AbstractMesh) = dofs.(nodes(m))

"Return true if the `AbstractMesh` m has `Dof`s defined."
_isempty_dofs(m::AbstractMesh) = !all(isempty.(dofs(m)))

"Return the number of `Dof`s defined in the `AbstractMesh` `m`.
This function assumes that `Dof`s indexes start from `Dof(1)`"
function num_dofs(m::AbstractMesh)::Int
    mesh_dofs = dofs(m)
    max_dof = 0
    !_isempty_dofs(m) && return max_dof # mesh has no dofs

    for node_dof in mesh_dofs
        for dofs in values(node_dof)
            max_dof = maximum([max_dof, maximum(dofs)])
        end
    end
    max_dof
end

"Adds n `dofs_per_node` `Dof`s with `dof_symbol` to the `AbstractMesh` `m`."
function apply!(m::AbstractMesh, dof_symbol::Symbol, dofs_per_node::Int)
    mesh_dofs = dofs(m)
    dof_not_added_yet = dof_symbol ∉ keys.(mesh_dofs)
    @assert dof_not_added_yet throw(ArgumentError("Dof symbol $dof_symbol already exists."))

    if !_isempty_dofs(m) && dof_not_added_yet  # any dof has been added
        for (i, n) in enumerate(nodes(m))
            node_dofs_int = (1 + (i - 1) * dofs_per_node):(i * dofs_per_node)
            apply!(n, dof_symbol, Dof.(node_dofs_int))
        end
    else # other dof has been added
        # Maximum dof index among all dofs
        max_dof_index = num_dofs(m)
        # Push new dofs
        for (i, n) in enumerate(nodes(m))
            node_dofs_int = (1 + max_dof_index + (i - 1) * dofs_per_node):(max_dof_index + i * dofs_per_node)
            apply!(n, dof_symbol, Dof.(node_dofs_int))
        end
    end
end

"Return a `Vector` of `Node`s defined in the `AbstractMesh` `m`."
nodes(m::AbstractMesh) = m.nodes

"Return the number of `Node`s of the `AbstractMesh` `m`."
num_nodes(m::AbstractMesh) = length(nodes(m))

"Return a `Vector` of `Face`s defined in the `AbstractMesh` `m`."
faces(m::AbstractMesh) = m.faces

"Return the number of `Face`s of the `AbstractMesh` `m`."
num_faces(m::AbstractMesh) = length(faces(m))

"Return the `Element`s of the `AbstractMesh` `m`."
elements(m::AbstractMesh) = m.elements

"Return the `Element`s of the `AbstractMesh` `m` at index `i`."
element(m::AbstractMesh, i::Int) = m.elements[i]

"Return the number of elements of the `AbstractMesh` `m`."
num_elements(m::AbstractMesh) = length(elements(m))

"Push a new `Node` into the `AbstractMesh` `m` and return the new `Node` position."
function Base.push!(m::AbstractMesh, n::AbstractNode)
    push!(nodes(m), n)
    length(nodes(m))
end

"Push a new  `Vector of `Node`s into the `AbstractMesh` `m`."
Base.push!(m::AbstractMesh, vn::Vector{<:AbstractNode}) = [push!(nodes(m), n) for n in vn]

"Push a new `Face` `f` into the `AbstractMesh` `m` and return the new `Face` position."
function Base.push!(m::AbstractMesh, f::AbstractFace)
    push!(faces(m), f)
    length(faces(m))
end

"Push a new `Vector` of `Face`s `vf` into the `AbstractMesh` `m`."
Base.push!(m::AbstractMesh, vf::Vector{<:AbstractFace}) = [push!(faces(m), f) for f in vf]

"Push a new `Element` into the `AbstractMesh` `m` and return the new `Element` position."
function Base.push!(m::AbstractMesh, e::AbstractElement)
    push!(elements(m), e)
    length(elements(m))
end

"Push a new vector of `Element`s into the `AbstractMesh` `m`."
Base.push!(m::AbstractMesh, ve::Vector{<:AbstractElement}) = [push!(elements(m), e) for e in ve]

""" 
A `Mesh` is a collection of `Element`s, `Face`s and `Node`s that cover the discretized domain, 
together with Sets of elements and nodes. 

### Methods:
* [`node_set`](@ref)
* [`add_node_to_set!`](@ref)
* [`element_set`](@ref)
* [`add_element_to_set!`](@ref)
* [`face_set`](@ref)
* [`add_face_to_set!`](@ref)
"""
struct Mesh{dim,N<:AbstractNode{dim},E<:AbstractElement,F<:AbstractFace,EX} <: AbstractMesh{dim}
    "`Node`s of the mesh with dimension `dim`."
    nodes::Vector{N}
    "`Element`s of the mesh."
    elements::Vector{E}
    "`Face`s of the mesh."
    faces::Vector{F}
    "Sets of `Node`s mapped to a string key."
    node_sets::Dictionary{String,Set{Int}}
    "Sets of `Element`s mapped to a string key."
    element_sets::Dictionary{String,Set{Int}}
    "Sets of `Face`s mapped to a string key."
    face_sets::Dictionary{String,Set{Int}}
    "Additional data or mesh info."
    extra::EX
    function Mesh(nodes::Vector{N},
                  elements::Vector{E}=Vector{AbstractElement}(),
                  faces::Vector{F}=Vector{AbstractFace}(),
                  node_sets=Dictionary{String,Set{Int}}(),
                  face_sets=Dictionary{String,Set{Int}}(),
                  element_sets=Dictionary{String,Set{Int}}(),
                  extra::EX=nothing) where
             {dim,N<:AbstractNode{dim},F<:AbstractFace,E<:AbstractElement,EX}
        new{dim,N,E,F,EX}(nodes, elements, faces, node_sets, element_sets, face_sets, extra)
    end
end

"Constructor for `Mesh` with a `Node`'s`Vector` and extra."
function Mesh(nodes::Vector{N}, extra::EX) where {dim,N<:AbstractNode{dim},EX}
    Mesh(nodes,
         Vector{AbstractElement}(),
         Vector{AbstractFace}(),
         Dictionary{String,Set{Int}}(),
         Dictionary{String,Set{Int}}(),
         Dictionary{String,Set{Int}}(),
         extra)
end

"Return the `Mesh` `m` `Node` `Set`s. "
node_set(m::Mesh) = m.node_sets

"Return the `Mesh` `m` `Node` `Set` with `node_set_name`. "
node_set(m::Mesh, node_set_name::S) where {S} = node_set(m)[String(node_set_name)]

"Return a `Vector` of `Node`s  with `node_set_name` in the `Mesh` `m`."
function nodes(m::Mesh, node_set_name::S) where {S}
    node_indexes = node_set(m, node_set_name)
    mesh_nodes = nodes(m)
    [mesh_nodes[i] for i in node_indexes]
end

"Add a `node_id` to the `Mesh` `m` `Set` `node_set_name`."
function add_node_to_set!(m::Mesh, node_set_name::S, node_id::Int) where {S}
    node_sets = node_set(m)
    if haskey(node_sets, node_set_name)
        push!(node_sets[String(node_set_name)], node_id)
    else
        insert!(node_sets, String(node_set_name), Set(node_id))
    end
    node_sets
end

"Add a `Node` `n` to the `Set` `node_set_name`."
function add_node_to_set!(m::Mesh, node_set_name::S, n::AbstractNode) where {S}
    node_id = findfirst(n_mesh -> n_mesh == n, nodes(m))
    isempty(node_id) && error("Node $n is not found in mesh.")
    add_node_to_set!(m, node_set_name, node_id)
end

"Return the `Mesh` `m` `Element` `Set`s. "
element_set(m::Mesh) = m.element_sets

"Return the `Mesh` `m` `Element` `Set` with `element_set_name`. "
element_set(m::Mesh, element_set_name::S) where {S} = element_set(m)[String(element_set_name)]

"Return a `Vector` of `Element`s  with `element_set_name` in the `Mesh` `m`."
function elements(m::Mesh, element_set_name::S) where {S}
    element_indexes = element_set(m, element_set_name)
    mesh_elements = elements(m)
    [mesh_elements[i] for i in element_indexes]
end

function Base.show(io::IO, m::Mesh)
    nnodes = length(m.nodes)
    nelems = length(m.elements)
    nfaces = length(m.faces)
    println("• Mesh with $nnodes nodes, $nelems elements and $nfaces faces.")
end

"Add `element_id` to the `Mesh` `m` `Set` `element_set_name`."
function add_element_to_set!(m::Mesh, element_set_name::S, element_id::Int) where {S}
    element_sets = element_set(m)
    if haskey(element_sets, element_set_name)
        push!(element_sets[String(element_set_name)], element_id)
    else
        insert!(element_sets, String(element_set_name), Set(element_id))
    end
    element_sets
end

"Add an `Element` `e` to  the `Set` `element_set_name`."
function add_element_to_set!(m::Mesh, element_set_name::S, e::AbstractElement) where {S}
    element_id = findfirst(e_mesh -> e_mesh == e, elements(m))
    isempty(element_id) && error("face $n is not found in mesh.")
    add_element_to_set!(m, element_set_name, element_id)
end

"Return the `Mesh` `m` `Face` `Set`s. "
face_set(m::Mesh) = m.face_sets

"Return the `Mesh` `m` `Element` `Set` with `element_set_name`. "
face_set(m::Mesh, fae_set_name::S) where {S} = face_set(m)[String(fae_set_name)]

"Return a `Vector` of `Face`s  with `face_set_name` in the `Mesh` `m`."
function faces(m::Mesh, face_set_name::S) where {S}
    face_indexes = face_set(m, face_set_name)
    mesh_faces = faces(m)
    [mesh_faces[i] for i in face_indexes]
end

"Add `face_id` to the `Mesh` `m` `Set` `face_set_name`."
function add_face_to_set!(m::Mesh, face_set_name::S, face_id::Int) where {S}
    face_sets = face_set(m)
    if haskey(face_sets, face_set_name)
        push!(face_sets[String(face_set_name)], face_id)
    else
        insert!(face_sets, String(face_set_name), Set(face_id))
    end
    face_sets
end

"Add a `Face` `n` to the `Set` `face_set_name`."
function add_face_to_set!(m::Mesh, face_set_name::S, f::AbstractFace) where {S}
    face_id = findfirst(f_mesh -> f_mesh == f, faces(m))
    isempty(face_id) && error("face $n is not found in mesh.")
    add_face_to_set!(m, face_set_name, face_id)
end

end # module
