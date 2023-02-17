export GenericCrossSection

""" Generic cross-section.

This generic cross-section struct can be used to define a cross-section not belonging to the existing categories. 

### Fields:
- `area` -- area.
- `Ixx` -- moment of inertia respect to `x` axis.
- `Iyy` -- moment of inertia respect to `y` axis.
- `Izz` -- moment of inertia respect to `z` axis.
- `Ixy` -- product moment of area respect to the `x-y` axes.
- `Ixz` -- product moment of area respect to the `x-z` axes.
- `Iyz` -- product moment of area respect to the `y-z` axes.
"""
struct GenericCrossSection{T<:Real} <: AbstractCrossSection
    area::T
    Ixx::T
    Iyy::T
    Izz::T
    Ixy::T
    Ixz::T
    Iyz::T
end

"Constructor of a `GenericCrossSection` if x, y and z are the principal axes."
GenericCrossSection(A, Ixx, Iyy, Izz) = GenericCrossSection(A, Ixx, Iyy, Izz, 0.0, 0.0, 0.0)
