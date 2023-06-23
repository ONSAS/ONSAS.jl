"Module defining frame elements."
module Frames

using Reexport
using StaticArrays

using ..CrossSections
using ..Entities
using ..HyperElasticMaterials
using ..Nodes
using ..Utils

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

function internal_forces(m::AbstractHyperElasticMaterial, f::Frame, u_e::AbstractVector)
    error("Not implemented.")
end

end
