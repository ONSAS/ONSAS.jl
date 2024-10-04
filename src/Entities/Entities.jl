"""
Module defining an interface, `AbstractEntity`.

Each entity consists of an object with material and nodes. The following entities are defined:

- Entities.
- Faces.
"""
module Entities

using Reexport, Dictionaries, StaticArrays

# using ..Nodes
using ..Materials

@reexport import ..Utils: label, apply!, dofs
@reexport import ..Nodes: coordinates, dimension
@reexport import ..CrossSections: area

export AbstractEntity, nodes, create_entity
export AbstractFace, normal_direction, volume
export AbstractElement, AbstractElementCache, cross_section, internal_forces, inertial_forces,
       local_dof_symbol, local_dofs, nodes, strain, stress, weights, num_nodes, elements_cache

# =================
# Abstract Entity
# =================
""" Abstract supertype for all `Face`s and `Element`s.

An `AbstractEntity` object is an entity defined by dofs and node/s with certain coordinates and dimension.

**Abstract Methods**

* [`create_entity`](@ref)
* [`coordinates`](@ref)
* [`dimension`](@ref)
* [`dofs`](@ref)

**Abstract fields**
* nodes
* label
"""
abstract type AbstractEntity{dim,T} end

"Return the `AbstractEntity` `e` coordinates."
coordinates(e::AbstractEntity) = coordinates.(nodes(e))

"Return each `AbstractFace` coordinates in a `Vector` of `Face`s `vf`."
coordinates(ve::Vector{<:AbstractEntity}) = coordinates.(ve)

"Return an `AbstractEntity` given an empty `AbstractEntity` `e` and a `Vector` of `Node`s `vn`."
# function create_entity(e::AbstractEntity, vn::AbstractVector{<:AbstractNode}) end

"Return the `AbstractEntity` dimension."
# dimension(::AbstractEntity{dim}) where {dim} = dim

"Return the dofs of an `AbstractEntity` `e`."
function dofs(e::AbstractEntity)
    mapfoldl(dofs, mergewith!(vcat), nodes(e); init=Dictionary{Symbol,Vector{Dof}}())
end

"Return the dofs of a `Vector` `ve` with `AbstractEntity`es."
dofs(ve::Vector{<:AbstractEntity}) = unique(row_vector(dofs.(ve)))

"Return the `Node`s an `AbstractEntity` `e`."
nodes(e::AbstractEntity) = e.nodes

"Return the label of an `AbstractEntity` `e`."
label(e::AbstractEntity) = e.label

# =================
# Abstract Face
# =================

""" Abstract supertype for all elements.

An `AbstractFace` object facilitates the process of adding boundary conditions on a surface.

**Abstract Methods**

* [`area`](@ref)
* [`create_entity`](@ref)
* [`coordinates`](@ref)
* [`dofs`](@ref)
* [`label`](@ref)
* [`nodes`](@ref)
* [`normal_direction`](@ref)

**Abstract fields**
* nodes
* label
"""
abstract type AbstractFace{dim,T} <: AbstractEntity{dim,T} end

"Return the `AbstractFace` `f` area."
function area(f::AbstractFace) end

"Return the `AbstractFace` `f` normal direction."
function normal_direction(f::AbstractFace) end

## =================
# Abstract Element
# =================

""" Abstract supertype for all elements.
An `AbstractElement` object facilitates the process of evaluating:
    - The internal forces vector and its tangent matrices.
    - The inertial forces vector and its tangent matrices.
    - The mechanical stresses and strains.

**Abstract Methods**
* [`coordinates`](@ref)
* [`dimension`](@ref)
* [`dofs`](@ref)
* [`local_dofs`](@ref)

These methods is a hard contract and must be implemented to define a new element.
* [`local_dof_symbol`](@ref)
* [`label`](@ref)
* [`nodes`](@ref)
* [`weights`](@ref)
* [`Base.:∈`](@ref)

This method is a hard contract and for static analysis must be implemented to define a new element.
* [`internal_forces`](@ref)

This method is a hard contract and for dynamic analysis must be implemented to define a new element.
* [`inertial_forces`](@ref)

**Abstract fields**
* nodes
* label
"""
abstract type AbstractElement{dim,T} <: AbstractEntity{dim,T} end

"Return true if a point `p` is inside the `AbstractElement` `e`."
function Base.:∈(p::AbstractVector, ::AbstractElement) end

"Return the `AbstractElement` `e` cross_section."
cross_section(e::AbstractElement) = e.cross_section

"Return local dofs symbols of the `AbstractElement` `e` (for linear displacements `:u` is used) in a vector.
Since global degrees of freedom are for the assemble process this function is used to compute the global dofs of the element by
extracting the node dofs with the symbol defined by the `AbstractElement` `e`."
function local_dof_symbol(e::AbstractElement) end

"""
Return local dofs given a vector of local dof symbols.
This method extracts all node dofs with the same symbol as local_dof_symbol.
"""
function local_dofs(e::AbstractElement)
    lds = local_dof_symbol(e)
    res = Dof[]
    for s in lds
        # Store in the resulting array the dofs per element that match each local dof.
        # Traversal order matters, since `local_dofs` is then used to build reduced matrices.
        for n in nodes(e)
            dict = dofs(n)
            if !haskey(dict, s)
                throw(ArgumentError("Element $(e.label) doesn't have dofs with symbol $s."))
            end
            append!(res, dict[s])
        end
    end
    res
end

"Return the internal forces vector of an `AbstractElement` `e` with an `AbstractMaterial` `m`."
function internal_forces(mat::AbstractMaterial, e::AbstractElement, u_e::AbstractVector)
    error("Not implemented.")
end

"Return the inertial forces vector of an `AbstractElement` `e`. with an `AbstractMaterial` `m`"
function inertial_forces(m::AbstractMaterial, e::AbstractElement, args...; kwargs...) end

"Return the `AbstractElement` `e` strain."
function strain(e::AbstractElement, args...; kwargs...) end

"Return the `AbstractElement` `e` stress."
function stress(e::AbstractElement, args...; kwargs...) end

"Return the weights to interpolate a scalar field at the `Node`s `Dof` corresponding
to the `AbstractElement` `e`."
function weights(e::AbstractElement, p::AbstractVector) end

"Return the number of nodes of element `e`."
function num_nodes(e::AbstractElement)
    length(nodes(e))
end

"Return the volume of the element."
function volume(e::AbstractElement) end

"""
Cache for internal forces computations.
"""
abstract type AbstractElementCache end

"Fallbaack cache for any element."
elements_cache(::Type{<:AbstractElement}) = nothing

"""
Concrete subtypes of `AbstractElement` should implement either the three-arg version (no cache)
or the four-arg version (use the cache).
"""
function internal_forces(mat::AbstractMaterial, e::AbstractElement, u_e::AbstractVector, ::Nothing)
    internal_forces(mat, e, u_e)
end
function internal_forces(mat::AbstractMaterial, e::AbstractElement, u_e::AbstractVector,
                         cache::AbstractElementCache)
    error("Not implemented.")
end

end # module
