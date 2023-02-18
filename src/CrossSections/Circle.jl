using ..CrossSections: AbstractCrossSection

import ..CrossSections: area, Ixx, Iyy, Izz, Ixy, Ixz, Iyz

export Circle

""" Circle cross-section.
### Fields:
- `d` -- circle diameter.
"""
struct Circle{T<:Real} <: AbstractCrossSection
    d::T
end

"Returns the area of a `Circle` cross-section `c`."
area(c::Circle) = pi * (c.d)^2 / 4

"Returns the moment of inertia of a `Circle` cross-section `c`  with respect to the local x axis."
Ixx(c::Circle) = pi * (c.d)^4 / 32

"Returns the moment of inertia of a `Circle` cross-section `c`  with respect to the local y axis."
Iyy(c::Circle) = pi * (c.d)^4 / 64

"Returns the moment of inertia of a `Circle` cross-section `c`  with respect to the local z axis."
Izz(c::Circle) = pi * (c.d)^4 / 64

"Returns the product moment of area of a `Circle` cross-section `c` with respect to the local x and y axes."
Ixy(c::Circle) = 0.0

"Returns the product moment of area of a `Circle` cross-section `c` with respect to the local x and z axes."
Ixz(c::Circle) = 0.0

"Returns the product moment of area of a `Circle` cross-section `c` with respect to the local y and z axes."
Iyz(c::Circle) = 0.0