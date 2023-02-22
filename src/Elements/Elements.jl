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

export Dof, add!
export AbstractNode, dimension, dofs, coordinates
export AbstractFace, coordinates, nodes
export AbstractElement, cross_section, internal_forces, inertial_forces, local_dof_symbol, local_dofs, nodes, strain, stress

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

"Returns the dof index from a vector of `Dof`s `vd` "
@inline index(vd::Vector{Dof}) = index.(vd)

"Returns the maximum dof index from a vector of `Dof`s `vd` "
Base.maximum(vd::Vector{Dof}) = maximum(index.(vd))

"Returns the entry of a vector `v` at index corresponding to the `Dof` `d` index."
@inline Base.getindex(v::AbstractVector, d::Dof) = v[index(d)]

"Sets the index of a vector `v` with the index corresponding to the `Dof` `d`."
@inline Base.setindex!(v::AbstractVector{T}, t::T, d::Dof) where {T} = v[index(d)] = t

"Returns a vector of entries of a vector `v` at indexes corresponding to the `Dof`s vector `vd`."
@inline Base.getindex(v::AbstractVector, vd::Vector{<:Dof}) = [v[index(d)] for d in vd]

"Sets a vector of values `tv` to a vector `v` at indexes corresponding to the `Dof`s vector `vd`."
@inline Base.setindex!(v::AbstractVector, tv::Vector{T}, vd::Vector{<:Dof}) where {T} = [setindex!(v, ti, vi) for (ti, vi) in zip(tv, vd)]

# =================
# Abstract Node
# =================

abstract type AbstractNode{dim,T} end

""" Abstract supertype for all nodes.

An `AbstractNode` object is a point in space.

**Common methods:**

* [`coordinates`](@ref)
* [`dimension`](@ref)
* [`dofs`](@ref)
"""

"Returns the `AbstractNode` `n` coordinates."
coordinates(n::AbstractNode) = n.x

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

"Sets a vectors of dofs `vd` to the `AbstractNode` `n` assigned to the symbol `s`."
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
# Abstract Face
# =================
abstract type AbstractFace{dim,T} end

""" Abstract supertype for all elements.

An `AbstractFace` object facilitates the process of adding boundary conditions on a surface. 

**Common methods:**

* [`coordinates`](@ref)
* [`dofs`](@ref)
* [`label`](@ref)
* [`nodes`](@ref)

**Common fields:**
* nodes
* label
"""

"Returns the `AbstractFace` `f` coordinates."
coordinates(f::AbstractFace) = coordinates.(nodes(f))

"Returns the `AbstractFace` `f` dimension."
dimension(::AbstractFace{dim}) where {dim} = dim

"Returns the dofs of `AbstractFace` `f`."
function dofs(f::AbstractFace)
    vecdfs = dofs.(nodes(f))
    dfs = mergewith(vcat, vecdfs[end], vecdfs[end-1])

    for i in 1:length(vecdfs)-2
        dfs = mergewith(vcat, dfs, vecdfs[end-(1+i)])
    end
    dfs
end

"Returns the dofs of a vector `vf` with `AbstractFace`s."
dofs(vf::Vector{<:AbstractFace}) = unique(row_vector(dofs.(vf)))

"Returns the label of `AbstractFace` `f`."
label(f::AbstractFace) = f.label

#==============================#
# AbstractFace implementations #
#==============================#

include("./TriangularFace.jl")


# =================
# Abstract Element
# =================

abstract type AbstractElement{dim,T} end

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
This method is a hard contract and must be implemented to define a new element.
* [`local_dof_symbol`](@ref)
* [`label`](@ref)
This method is a hard contract and must be implemented to define a new element.
* [`nodes`](@ref)

This method is a hard contract and for static analysis must be implemented to define a new element.
* [`internal_forces`](@ref)
This method is a hard contract and for dynamic analysis must be implemented to define a new element.
* [`inertial_forces`](@ref)

**Common fields:**
* nodes
* label

"""

"Returns the `AbstractElement` `e` coordinates."
coordinates(e::AbstractElement) = coordinates.(nodes(e))

"Returns the `AbstractElement` `e` cross_section."
cross_section(e::AbstractElement) = e.cross_section

"Returns the `AbstractElement` `e` dimension."
dimension(::AbstractElement{dim}) where {dim} = dim

"Returns the dofs of `AbstractElement` `e`."
dofs(e::AbstractElement) = mergewith(vcat, dofs.(nodes(e))...)

"Returns the dofs of a vector `ve` with `AbstractElement`s."
dofs(ve::Vector{<:AbstractElement}) = unique(row_vector(dofs.(ve)))

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


"Returns the label of an `AbstractElement` `e`."
label(e::AbstractElement) = e.label

"Returns the `Node`s of an `AbstractElement` `e`."
nodes(e::AbstractElement) = e.nodes

"Returns the internal forces vector of an `AbstractElement` `e`."
function internal_forces(e::AbstractElement, args...; kwargs...) end

"Returns the inertial forces vector of an `AbstractElement` `e`."
function inertial_forces(e::AbstractElement, args...; kwargs...) end

"Returns the `AbstractElement` `e` strain."
function strain(e::AbstractElement, args...; kwargs...) end

"Returns the `AbstractElement` `e` stress."
function stress(e::AbstractElement, args...; kwargs...) end

#=================================#
# AbstractElement implementations #
#=================================#

include("./Truss.jl")

end # module


