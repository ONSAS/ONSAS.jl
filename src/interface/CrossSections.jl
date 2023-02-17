"""
Module defining geometric entities interface.
Each cross section consists of a data type with one or more geometric parameters into its fields.
"""
module CrossSections

export AbstractCrossSection, area, Ixx, Iyy, Izz, Ixy, Ixz, Iyz, I_tensor_xyz
export Rectangle, Square, Circle, GenericCrossSection

# ======================
# Abstract Cross Section
# ======================

""" Abstract supertype for all elements corss-section.

**Common methods:**
The following functions must be implemented for a new cross-section:
- * [`area`](@ref)         -- returns the cross-section area.
- * [`Ixx`](@ref)          -- returns the moment of inertia with respect to `x` axis.
- * [`Iyy`](@ref)          -- returns the moment of inertia with respect to `y` axis.
- * [`Izz`](@ref)          -- returns the moment of inertia with respect to `z` axis.
- * [`Ixy`](@ref)          -- returns the product moment of area respect to the `x-y` axes.
- * [`Ixz`](@ref)          -- returns the product moment of area respect to the `x-z` axes.
- * [`Iyz`](@ref)          -- returns the product moment of area respect to the `y-z` axes.
- * [`I_tensor_xyz`](@ref) -- returns the inertia tensor in the system of coordinates `x-y-z`.
"""

abstract type AbstractCrossSection end

"Returns the cross-section area of the `AbstractCrossSection` `cs`"
area(cs::AbstractCrossSection) = cs.area

"Returns the moment of inertia of the `AbstractCrossSection` `cs` with respect to the local x axis."
Ixx(::AbstractCrossSection) = cs.Ixx

"Returns the moment of inertia of the `AbstractCrossSection` `cs` with respect to the local `y` axis."
Iyy(::AbstractCrossSection) = cs.Iyy

"Returns the moment of inertia of the `AbstractCrossSection` `cs` respect  to the local `z` axis."
Izz(::AbstractCrossSection) = cs.Izz

"Returns the product moment of area  of the `AbstractCrossSection` `cs` respect to the local `x-y` axes."
Ixy(::AbstractCrossSection) = cs.Ixy

"Returns the product moment of area of the `AbstractCrossSection` `cs` respect to the `x-z` local axes."
Ixz(::AbstractCrossSection) = cs.Ixz

"Returns the product moment of area of the `AbstractCrossSection` `cs` respect to the `y-z` axes."
Iyz(::AbstractCrossSection) = cs.Iyz

"Returns the inertia tensor in the system of coordinates `x-y-z`"
I_tensor_xyz(cs::AbstractCrossSection) = [
    Ixx(cs) -Ixy(cs) -Ixz(cs)
    -Ixy(cs) Iyy(cs) -Iyz(cs)
    -Ixz(cs) -Iyz(cs) Izz(cs)
]

#======================================#
# AbstractCrossSection implementations #
#======================================#

include("../cross_sections/Circle.jl")
include("../cross_sections/Rectangle.jl")
include("../cross_sections/Square.jl")
include("../cross_sections/GenericCrossSection.jl")

end # module





