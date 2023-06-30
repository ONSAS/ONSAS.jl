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
    # [u1_1, u2_1, u3_1, u1_2, u2_2, u3_2, t1_1, t2_1, t3_1, t1_2, t2_2, t3_2]

    # Temporalmente cero.
    σ = 0.0
    ε = 0.0

    E = elasticity_modulus(m)
    nu = poisson_ratio(m)
    G = shear_modulus(m)
    S = cross_section(f)
    A = area(S)
    J = CrossSections.Ixx(S)
    Iyy = CrossSections.Iyy(S)
    @show Izz = CrossSections.Izz(S)

    l = norm(f.nodes[2] - f.nodes[1])

    (ux1, uy1, uz1, ux2, uy2, uz2, titax1, titay1, titaz1, titax2, titay2, titaz2) = u_e

    Kloc = E * Izz / l^3 * [  12     6*l    -12     6*l
                              6*l   4*l^2   -6*l   2*l^2
                            -12    -6*l     12    -6*l
                              6*l   2*l^2   -6*l   4*l^2]

    Ks = zeros(12, 12)
    fint = zeros(12)

    # Bending along x-y.
    ind = [2, 9, 5, 12]
    Ks[ind, ind] .= Kloc
    fint .= Ks * u_e

    # fxx = Kloc * [uy1, titaz1, uy2, titaz2]
    return fint, Ks, σ, ε
end

end # module
