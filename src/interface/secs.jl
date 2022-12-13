abstract type AbstractSection end

struct Rectangle <: AbstractSection
    width_y::Float64
    width_z::Float64
end

struct Square <: AbstractSection
    dim::Float64
end

mutable struct CrossSection
    section::AbstractSection
    A::Float64
    Ix::Float64
    Iy::Float64
    Iz::Float64
end

function CrossSection(section::Rectangle)
    A = section.width_y * section.width_z
    Ix = 1
    Iy = 1
    Iz = 1
    return CrossSection(section, A, Ix, Iy, Iz)
end

function CrossSection(section::Square)
    A = section.dim
    Ix = 1
    Iy = 1
    Iz = 1
    return CrossSection(section, A, Ix, Iy, Iz)
end


