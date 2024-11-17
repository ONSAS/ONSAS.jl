module GenericCrossSections

using ..CrossSections: AbstractCrossSection

export GenericCrossSection

"""
This generic cross-section struct can be used to define a cross-section not belonging to the existing categories. 
"""
struct GenericCrossSection{T <: Real} <: AbstractCrossSection
    "Cross section area."
    area::T
    "Moment of inertia respect to `x` axis."
    Ixx::T
    "Moment of inertia respect to `y` axis."
    Iyy::T
    "Moment of inertia respect to `z` axis."
    Izz::T
    "Product moment of area respect to the `x-y` axes."
    Ixy::T
    "Product moment of area respect to the `x-z` axes."
    Ixz::T
    "Product moment of area respect to the `y-z` axes."
    Iyz::T
end

"Constructor of a `GenericCrossSection` if x, y and z are the principal axes."
GenericCrossSection(A, Ixx, Iyy, Izz) = GenericCrossSection(A, Ixx, Iyy, Izz, 0.0, 0.0, 0.0)

end # module
