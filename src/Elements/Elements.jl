"""
Module defining element formulations implemented.
Each element consists of a data type with the geometric properties such as nodes, cross-sections and a label into its fields.
"""
module Elements

using AutoHashEquals: @auto_hash_equals
using Reexport: @reexport
@reexport using Dictionaries: Dictionary, dictionary
using StaticArrays: Size, SVector, StaticArray
using ..Utils: row_vector

@reexport using ..Materials
@reexport using ..CrossSections
@reexport import ..Utils: label
@reexport import Dictionaries: index

import StaticArrays
import ..CrossSections: area

export Dof, add!, Point
export AbstractNode, dimension, dofs, coordinates
export AbstractEntity, nodes, coordinates, create_entity
export AbstractFace, normal_direction
export AbstractElement, cross_section, internal_forces, inertial_forces, local_dof_symbol, local_dofs, nodes, strain, stress, weights

const Point{dim,T} = Union{<:AbstractVector{P},NTuple{dim,P}} where {dim,P<:Real}

# ========================
# Degree of freedom (Dof)
# ========================

"""
Scalar degree of freedom of the structure.
"""
const Dof = Int

"Return the dof index of the `Dof` `d` "
index(d::Dof) = d

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
abstract type AbstractNode{dim,T} <: StaticArray{Tuple{dim},T,1} end

"Return the `AbstractNode` `n` coordinates."
coordinates(n::AbstractNode) = n.x

"Return each `AbstractNode` coordinates in a `Vector` of `Node`s vn."
coordinates(vn::AbstractVector{<:AbstractNode}) = coordinates.(vn)

"Return the `AbstractNode` `n` dimension (1D, 2D or 3D)."
dimension(::AbstractNode{dim}) where {dim} = dim

"Return a tuple containing the dimensions of the `AbstractNode` `n`."
Base.size(n::AbstractNode) = size(n.x)

"Return the coordinate at component `i` from the `AbstractNode` `n`."
Base.getindex(n::AbstractNode, i::Int) = n.x[i]

"Return `AbstractNode` `n` degrees of freedom."
dofs(n::AbstractNode) = n.dofs

"Return degrees of freedom for each `AbstractNode` in a vector of nodes `vn`."
dofs(vn::Vector{<:AbstractNode}) = vcat(dofs.(vn)...)

"Return `AbstractNode` `n` degrees of freedom with symbol `s`."
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

"Return the `AbstractEntity` `e` coordinates."
coordinates(e::AbstractEntity) = coordinates.(nodes(e))

"Return each `AbstractFace` coordinates in a `Vector` of `Face`s `vf`."
coordinates(ve::Vector{<:AbstractEntity}) = coordinates.(ve)

"Return an `AbstractEntity` given an empty `AbstractEntity` `e` and a `Vector` of `Node`s `vn`."
function create_entity(e::AbstractEntity, vn::AbstractVector{<:AbstractNode}) end

"Return the `AbstractEntity` dimension."
dimension(::AbstractEntity{dim}) where {dim} = dim

"Return the dofs of an `AbstractEntity` `e`."
function dofs(e::AbstractEntity)
    vecdfs = dofs.(nodes(e))
    dfs = mergewith(vcat, vecdfs[1], vecdfs[2])
    [mergewith!(vcat, dfs, vecdfs[i]) for i in 3:length(vecdfs)]
    dfs
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

"Return the `AbstractFace` `f` area."
function area(f::AbstractFace) end

"Return the `AbstractFace` `f` normal."
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
* [`Base.:∈`](@ref)

This method is a hard contract and for static analysis must be implemented to define a new element.
* [`internal_forces`](@ref)

This method is a hard contract and for dynamic analysis must be implemented to define a new element.
* [`inertial_forces`](@ref)

**Common fields:**
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

"Return local dofs given a vector of local dof symobls. This method extracts all node dofs with the same symbol
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

"Return the internal forces vector of an `AbstractElement` `e` with an `AbstractMaterial` `m`."
function internal_forces(m::AbstractMaterial, e::AbstractElement, args...; kwargs...) end

"Return the inertial forces vector of an `AbstractElement` `e`. with an `AbstractMaterial` `m`"
function inertial_forces(m::AbstractMaterial, e::AbstractElement, args...; kwargs...) end

"Return the `AbstractElement` `e` strain."
function strain(e::AbstractElement, args...; kwargs...) end

"Return the `AbstractElement` `e` stress."
function stress(e::AbstractElement, args...; kwargs...) end

"Return the weights to interpolate a scalar field at the `Node`s `Dof` corresponding 
to the `AbstractElement` `e`."
function weights(e::AbstractElement, p::AbstractVector) end

#=================================#
# AbstractElement implementations #
#=================================#

include("./Truss.jl")
include("./Tetrahedron.jl")

end # module


