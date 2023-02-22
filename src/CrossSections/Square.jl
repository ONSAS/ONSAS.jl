using ..CrossSections: AbstractCrossSection

import ..CrossSections: area, Ixx, Iyy, Izz, Ixy, Ixz, Iyz

export Square

""" Square cross-section.
### Fields:
- `width` -- width in `y` and `z` axes.
"""
struct Square{T<:Real} <: AbstractCrossSection
    width::T
end

area(s::Square) = s.width^2

"Returns the moment of inertia of a `Square` cross-section `s` with respect to the local x axis."
Ixx(s::Square) = (0.5 * s.width)^4 * (16 / 3 - 3.36 * (1 - (0.5 * s.width)^4 / (12 * (0.5 * s.width)^4)))

"Returns the moment of inertia of a `Square` cross-section `s` with respect to the local y axis."
Iyy(s::Square) = s.width^4 / 12

"Returns the moment of inertia of a `Square` cross-section `s` with respect to the local z axis."
Izz(s::Square) = s.width^4 / 12

"Returns the product moment of area of a `Square` cross-section `s` with respect to the local x and y axes."
Ixy(s::Square) = 0.0

"Returns the product moment of area of a `Square` cross-section `s` with respect to the local x and z axes."
Ixz(s::Square) = 0.0

"Returns the product moment of area of a `Square` cross-section `s` with respect to the local y and z axes."
Iyz(s::Square) = 0.0