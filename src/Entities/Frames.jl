"Module defining frame elements."
module Frames

using Reexport
using StaticArrays
using LinearAlgebra

using ..CrossSections
using ..Entities
using ..Nodes
using ..Utils
using ..IsotropicLinearElasticMaterial

@reexport import ..Entities: nodes, cross_section, create_entity, local_dof_symbol, internal_forces

export Frame, Consistent, Lumped

"Type of mass matrix."
@enum MassMatrix Consistent Lumped

struct Frame{dim,T<:Real,N<:AbstractNode{dim,T},VN<:AbstractVector{N},G<:AbstractCrossSection} <:
       AbstractElement{dim,T}
    "Nodes of the frame."
    nodes::VN
    "Cross section properties."
    cross_section::G
    "Type of mass matrix."
    mass_matrix::MassMatrix
    "Label of the frame."
    label::Label
    function Frame(nodes::VN, cross_section::G, mass_matrix::MassMatrix=Consistent,
                   label::Label=NO_LABEL) where {dim,T<:Real,N<:AbstractNode{dim,T},
                                                 VN<:AbstractVector{N},G<:AbstractCrossSection}
        @assert length(nodes) == 2 ||
                throw(ArgumentError("Expected two nodes, got $(length(nodes))"))
        new{dim,T,N,VN,G}(nodes, cross_section, mass_matrix, Symbol(label))
    end
end
function Frame(n₁::N, n₂::N, g::G, mass_matrix::MassMatrix=Consistent,
               label::Label=NO_LABEL) where {N<:AbstractNode,G<:AbstractCrossSection}
    Frame(SVector(n₁, n₂), g, mass_matrix, label)
end

nodes(f::Frame) = f.nodes
cross_section(f::Frame) = f.cross_section

function create_entity(f::Frame, vn::AbstractVector{<:AbstractNode})
    Frame(vn, cross_section(t), mass_matrix(t), label(t))
end

local_dof_symbol(::Frame) = [:u, :θ]

function internal_forces(m::IsotropicLinearElastic, f::Frame, u_e::AbstractVector)
    @assert f.mass_matrix == Consistent

    # Reference:
    # (ux1, uy1, uz1, ux2, uy2, uz2, θx1, θy1, θz1, θx2, θy2, θz2) = u_e

    # TODO Generalize.
    σ = 0.0
    ε = 0.0

    E = elasticity_modulus(m)
    G = shear_modulus(m)
    S = cross_section(f)
    A = area(S)
    J = CrossSections.Ixx(S)
    Iyy = CrossSections.Iyy(S)
    Izz = CrossSections.Izz(S)
    l = norm(f.nodes[2] - f.nodes[1])

    ind_bend_xy = [2, 9, 5, 12]
    ind_bend_xz = [3, 8, 6, 11]
    inds_axial = [1, 4]
    inds_torsion = [7, 10]

    Ks = zeros(12, 12)
    fint = zeros(12)

    Kbend = [12     6*l    -12     6*l
             6*l   4*l^2   -6*l   2*l^2
             -12    -6*l     12    -6*l
             6*l   2*l^2   -6*l   4*l^2]

    # Bending along x-y.
    Ks[ind_bend_xy, ind_bend_xy] .+= E * Izz / l^3 * Kbend

    # Bending along x-z.
    Rxyxz = diagm([1, -1, 1, -1])
    Ks[ind_bend_xz, ind_bend_xz] .+= E * Iyy / l^3 * Rxyxz * Kbend * Rxyxz

    # Axial stiffness along x.
    Ks[inds_axial, inds_axial] .+= E * A / l * [1 -1; -1 1]

    # Torsion stiffness along x.
    Ks[inds_axial, inds_torsion] .+= G * J / l * [1 -1; -1 1]

    # Internal forces.
    fint .= Ks * u_e

    return fint, Ks, σ, ε
end

end # module
