"""
Module defining geometric entities interface.
Each cross section consists of a data type with one or more geometric parameters into its fields.
"""
module CrossSections

export AbstractCrossSection, area, Ixx, Iyy, Izz, Ixy, Ixz, Iyz, I_tensor_xyz

# ======================
# Abstract Cross Section
# ======================

"""
Abstract supertype for all elements corss-section.

## Common methods

The following functions must be implemented for a new cross-section:
- * [`area`](@ref)         -- Return the cross-section area.
- * [`Ixx`](@ref)          -- Return the moment of inertia with respect to `x` axis.
- * [`Iyy`](@ref)          -- Return the moment of inertia with respect to `y` axis.
- * [`Izz`](@ref)          -- Return the moment of inertia with respect to `z` axis.
- * [`Ixy`](@ref)          -- Return the product moment of area respect to the `x-y` axes.
- * [`Ixz`](@ref)          -- Return the product moment of area respect to the `x-z` axes.
- * [`Iyz`](@ref)          -- Return the product moment of area respect to the `y-z` axes.
- * [`I_tensor_xyz`](@ref) -- Return the inertia tensor in the system of coordinates `x-y-z`.
"""
abstract type AbstractCrossSection end

"Return the cross-section area of the `AbstractCrossSection` `cs`"
area(cs::AbstractCrossSection) = cs.area

"Return the moment of inertia of the `AbstractCrossSection` `cs` with respect to the local x axis."
Ixx(::AbstractCrossSection) = cs.Ixx

"Return the moment of inertia of the `AbstractCrossSection` `cs` with respect to the local `y` axis."
Iyy(::AbstractCrossSection) = cs.Iyy

"Return the moment of inertia of the `AbstractCrossSection` `cs` respect  to the local `z` axis."
Izz(::AbstractCrossSection) = cs.Izz

"Return the product moment of area  of the `AbstractCrossSection` `cs` respect to the local `x-y` axes."
Ixy(::AbstractCrossSection) = cs.Ixy

"Return the product moment of area of the `AbstractCrossSection` `cs` respect to the `x-z` local axes."
Ixz(::AbstractCrossSection) = cs.Ixz

"Return the product moment of area of the `AbstractCrossSection` `cs` respect to the `y-z` axes."
Iyz(::AbstractCrossSection) = cs.Iyz

"Return the inertia tensor in the system of coordinates `x-y-z`"
function I_tensor_xyz(cs::AbstractCrossSection)
    return [Ixx(cs) -Ixy(cs) -Ixz(cs)
            -Ixy(cs) Iyy(cs) -Iyz(cs)
            -Ixz(cs) -Iyz(cs) Izz(cs)]
end

end # module
