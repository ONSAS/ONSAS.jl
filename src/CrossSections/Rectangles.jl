module Rectangles

using ..CrossSections: AbstractCrossSection

import ..CrossSections: area, Ixx, Iyy, Izz, Ixy, Ixz, Iyz

export Rectangle

"""
Rectangle cross-section.
"""
struct Rectangle{T<:Real} <: AbstractCrossSection
    "Width in `y` local axis."
    width_y::T
    "Width in `z` local axis."
    width_z::T
end

"Return the area of a `Rectangle` cross-section `r`."
area(r::Rectangle) = r.width_y * r.width_z

"Return the moment of inertia of a `Rectangle` cross-section `r`  with respect to the local x axis."
function Ixx(r::Rectangle)
    #  torsional constant from table 10.1 from Roark's Formulas for Stress and Strain 7th ed.
    a = 0.5 * maximum([r.width_y, r.width_z])
    b = 0.5 * minimum([r.width_y, r.width_z])

    Ixx = a * b^3 * (16 / 3 - 3.36 * b / a * (1 - b^4 / (12 * a^4)))
    return Ixx
end

"Return the moment of inertia of a `Rectangle` cross-section `r`  with respect to the local y axis."
Iyy(r::Rectangle) = r.width_z^3 * r.width_y / 12

"Return the moment of inertia of a `Rectangle` cross-section `r`  with respect to the local z axis."
Izz(r::Rectangle) = r.width_y^3 * r.width_z / 12

"Return the product moment of area of a `Rectangle` cross-section `r` with respect to the local x and y axes."
Ixy(r::Rectangle) = 0.0

"Return the product moment of area of a `Rectangle` cross-section `r` with respect to the local x and z axes."
Ixz(r::Rectangle) = 0.0

"Return the product moment of area of a `Rectangle` cross-section `r` with respect to the local y and z axes."
Iyz(r::Rectangle) = 0.0

end # module
