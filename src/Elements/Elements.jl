"""
Module defining element formulations implemented.
Each element consists of a data type with the geometric properties such as nodes, cross-sections and a label into its fields.
"""
module Elements

using AutoHashEquals: @auto_hash_equals
using Reexport: @reexport
@reexport using Dictionaries: Dictionary, dictionary
using StaticArrays: SVector
using ..Utils: row_vector

@reexport using ..Materials
@reexport using ..CrossSections
@reexport import ..Utils: label
@reexport import Dictionaries: index

import ..CrossSections: area

export Dof, add!
export AbstractNode, dimension, dofs, coordinates
export AbstractEntity, nodes, coordinates, create_entity
export AbstractFace, normal_direction
export AbstractElement, cross_section, internal_forces, inertial_forces, local_dof_symbol, local_dofs, nodes, strain, stress, weights

# ========================
# Degree of freedom (Dof)
# ========================

""" Degree of freedom struct.
This is a scalar degree of freedom of the structure.
### Fields:
- `index`   -- degree of freedom identification number. 
"""
@auto_hash_equals struct Dof
    index::Int
end

"Returns the dof index of the `Dof` `d` "
@inline index(d::Dof) = d.index

"Returns the dof index from a `Vector` of `Dof`s `vd` "
@inline index(vd::Vector{Dof}) = index.(vd)

"Returns the maximum dof index from a `Vector` of `Dof`s `vd` "
Base.maximum(vd::Vector{Dof}) = maximum(index.(vd))

"Returns the entry of a `Vector` `v` at index corresponding to the `Dof` `d` index."
@inline Base.getindex(v::AbstractVector, d::Dof) = v[index(d)]

"Sets the index of a `Vector` `v` with the index corresponding to the `Dof` `d`."
@inline Base.setindex!(v::AbstractVector{T}, t::T, d::Dof) where {T} = v[index(d)] = t

"Returns a `Vector` of entries of a `Vector` `v` at indexes corresponding to the `Dof`s vector `vd`."
@inline Base.getindex(v::AbstractVector, vd::Vector{<:Dof}) = [v[index(d)] for d in vd]

"Sets a `Vector` of values `tv` to a `Vector` `v` at indexes corresponding to the `Dof`s vector `vd`."
@inline Base.setindex!(v::AbstractVector, tv::Vector{T}, vd::Vector{<:Dof}) where {T} = [setindex!(v, ti, vi) for (ti, vi) in zip(tv, vd)]

# =================
# Abstract Node
# =================

""" Abstract supertype for all nodes.

An `AbstractNode` object is a point in space.

**Common methods:**

* [`coordinates`](@ref)
* [`dimension`](@ref)
* [`dofs`](@ref)
"""
abstract type AbstractNode{dim,T} end

"Returns the `AbstractNode` `n` coordinates."
coordinates(n::AbstractNode) = n.x

"Returns each `AbstractNode` coordinates in a `Vector` of `Node`s vn."
coordinates(vn::Vector{<:AbstractNode}) = coordinates.(vn)

"Returns the `AbstractNode` `n` dimension (1D, 2D or 3D)."
dimension(::AbstractNode{dim}) where {dim} = dim

"Returns the coordinate at component `i` from the `AbstractNode` `n`."
Base.getindex(n::AbstractNode, i::Int) = n.x[i]

"Returns `AbstractNode` `n` degrees of freedom."
dofs(n::AbstractNode) = n.dofs

"Returns degrees of freedom for each `AbstractNode` in a vector of nodes `vn`."
dofs(vn::Vector{<:AbstractNode}) = vcat(dofs.(vn)...)

"Returns `AbstractNode` `n` degrees of freedom with symbol `s`."
dofs(n::AbstractNode, s::Symbol) = n.dofs[s]

"Sets a `Vector`s of dofs `vd` to the `AbstractNode` `n` assigned to the symbol `s`."
function add!(n::AbstractNode, s::Symbol, vd::Vector{Dof})
    if s ∉ keys(dofs(n))
        insert!(dofs(n), s, vd)
    else
        [push!(dofs(n)[s], d) for d in vd if d ∉ dofs(n)[s]]
    end
    return dofs(n)
end

include("./Node.jl")

# =================
# Abstract Entity
# =================
""" Abstract supertype for all `Face`s and `Element`s.

An `AbstractEntity` object is an entity defined by dofs and node/s with certain coordinates and dimension.

**Common methods:**

* [`create_entity`](@ref)
* [`coordinates`](@ref)
* [`dimension`](@ref)
* [`dofs`](@ref)

**Common fields:**
* nodes
* label
"""
abstract type AbstractEntity{dim,T} end

"Returns the `AbstractEntity` `e` coordinates."
coordinates(e::AbstractEntity) = coordinates.(nodes(e))

"Returns each `AbstractFace` coordinates in a `Vector` of `Face`s `vf`."
coordinates(ve::Vector{<:AbstractEntity}) = coordinates.(ve)

"Returns an `AbstractEntity` given an empty `AbstractEntity` `e` and a `Vector` of `Node`s `vn`."
function create_entity(e::AbstractEntity, vn::AbstractVector{<:AbstractNode}) end

"Returns the `AbstractEntity` dimension."
dimension(::AbstractEntity{dim}) where {dim} = dim

"Returns the dofs of an `AbstractEntity` `e`."
function dofs(e::AbstractEntity)
    vecdfs = dofs.(nodes(e))
    dfs = mergewith(vcat, vecdfs[1], vecdfs[2])
    [mergewith!(vcat, dfs, vecdfs[i]) for i in 3:length(vecdfs)]
    dfs
end

"Returns the dofs of a `Vector` `ve` with `AbstractEntity`es."
dofs(ve::Vector{<:AbstractEntity}) = unique(row_vector(dofs.(ve)))

"Returns the `Node`s an `AbstractEntity` `e`."
nodes(e::AbstractEntity) = e.nodes

"Returns the label of an `AbstractEntity` `e`."
label(e::AbstractEntity) = e.label

# =================
# Abstract Face
# =================

""" Abstract supertype for all elements.

An `AbstractFace` object facilitates the process of adding boundary conditions on a surface. 

**Common methods:**

* [`area`](@ref)
* [`create_entity`](@ref)
* [`coordinates`](@ref)
* [`dofs`](@ref)
* [`label`](@ref)
* [`nodes`](@ref)
* [`normal_direction`](@ref)

**Common fields:**
* nodes
* label
"""
abstract type AbstractFace{dim,T} <: AbstractEntity{dim,T} end

"Returns the `AbstractFace` `f` area."
function area(f::AbstractFace) end

"Returns the `AbstractFace` `f` normal."
function normal_direction(f::AbstractFace) end

#==============================#
# AbstractFace implementations #
#==============================#

include("./TriangularFace.jl")

## =================
# Abstract Element
# =================

""" Abstract supertype for all elements.
An `AbstractElement` object facilitates the process of evaluating:
    - The internal forces vector and its tangent matrices.
    - The inertial forces vector and its tangent matrices.
    - The mechanical stresses and strains.

**Common methods:**
* [`coordinates`](@ref)
* [`dimension`](@ref)
* [`dofs`](@ref)
* [`local_dofs`](@ref) 

These methods is a hard contract and must be implemented to define a new element.
* [`local_dof_symbol`](@ref)
* [`label`](@ref)
* [`nodes`](@ref)
* [`weights`](@ref)

This method is a hard contract and for static analysis must be implemented to define a new element.
* [`internal_forces`](@ref)

This method is a hard contract and for dynamic analysis must be implemented to define a new element.
* [`inertial_forces`](@ref)

**Common fields:**
* nodes
* label
"""
abstract type AbstractElement{dim,T} <: AbstractEntity{dim,T} end

"Returns the `AbstractElement` `e` cross_section."
cross_section(e::AbstractElement) = e.cross_section

"Returns local dofs symbols of the `AbstractElement` `e` (for linear displacements `:u` is used) in a vector.
Since global degrees of freedom are for the assemble process this function is used to compute the global dofs of the element by 
extracting the node dofs with the symbol defined by the `AbstractElement` `e`."
function local_dof_symbol(e::AbstractElement) end

"Returns local dofs given a vector of local dof symobls. This method extracts all node dofs with the same symbol
as local_dof_symbol"
function local_dofs(e::AbstractElement)
    local_dof_symbols = local_dof_symbol(e)
    local_dofs = Vector{Dof}()
    element_dofs = dofs(e)
    for dof_symbol in local_dof_symbols
        if dof_symbol ∉ keys(element_dofs)
            error("Element $(e.label) does not have dofs with symbol $(dof_symbol)")
        else
            push!(local_dofs, element_dofs[dof_symbol]...)
        end
    end
    return local_dofs
end

"Returns the internal forces vector of an `AbstractElement` `e` with an `AbstractMaterial` `m`."
function internal_forces(m::AbstractMaterial, e::AbstractElement, args...; kwargs...) end

"Returns the inertial forces vector of an `AbstractElement` `e`. with an `AbstractMaterial` `m`"
function inertial_forces(m::AbstractMaterial, e::AbstractElement, args...; kwargs...) end

"Returns the `AbstractElement` `e` strain."
function strain(e::AbstractElement, args...; kwargs...) end

"Returns the `AbstractElement` `e` stress."
function stress(e::AbstractElement, args...; kwargs...) end

"Returns the weights to interpolate a scalar field at the `Node`s `Dof` corresponding 
to the `AbstractElement` `e`."
function weights(e::AbstractElement, p::Any, args...; kwargs...) end

#=================================#
# AbstractElement implementations #
#=================================#

include("./Truss.jl")
include("./Tetrahedron.jl")

end # module


