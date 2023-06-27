"Module defining frame elements."
module Frames

using Reexport
using StaticArrays

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
    @show m
    @show f
    @show u_e

    # material constit params
    E  = elasticity_modulus(m)
    nu = poisson_ratio(m)
    G  = shear_modulus(m)

    S = cross_section(f)
    A = area(S)
    J = CrossSections.Ixx(S)
    Iyy = CrossSections.Iyy(S)
    Izz = CrossSections.Izz(S)

    #     # --- elem lengths and rotation matrix
    # 	[ local2globalMats, l ] = beamParameters( elemCoords ) ;
    # 	R = RotationMatrix(ndofpnode, local2globalMats) ;

    #   % temporary
    #   %~ ------------------------
    #   elemReleases = [0 0 0 0] ;
    #   %~ ------------------------

    #   % --- set the local degrees of freedom corresponding to each behavior
    #   LocAxialdofs  = [ 1 7 ] ;
    #   LocTorsndofs  = [ 2 8 ] ;
    #   LocBendXYdofs = [ 3 6 9 12 ] ;
    #   LocBendXZdofs = [ 5 4 11 10 ] ;

    #   KL = zeros ( 2*ndofpnode, 2*ndofpnode ) ;

    #   Kaxial = E*A/l * [ 1 -1 ; ...
    #                     -1  1 ] ;
    #   KL( LocAxialdofs , LocAxialdofs ) = Kaxial ;

    #   kBendReleaseRig = [ 3    3*l   -3   0 ; ...
    #                       3*l  3*l^2 -3*l 0 ; ...
    #                      -3   -3*l    3   0 ; ...
    #                       0    0      0   0 ] ;

    #   kBendReleaseLef = [  3   0 -3   3*l   ; ...
    #                        0   0  0   0     ; ...
    #                       -3   0  3  -3*l   ; ...
    #                        3*l 0 -3*l 3*l^2 ] ;

    #   % K bending in local coordinates
    #   kBendNoRelease = [  12     6*l    -12     6*l   ; ...
    #                        6*l   4*l^2   -6*l   2*l^2 ; ...
    #                      -12    -6*l     12    -6*l   ; ...
    #                        6*l   2*l^2   -6*l   4*l^2 ] ;

    #   % bending XY
    #   if     elemReleases(3) == 0 && elemReleases(4) == 0
    #     KbendXY = E * Iz / l^3 * kBendNoRelease ;
    #   elseif elemReleases(3) == 1 && elemReleases(4) == 0
    #     KbendXY = E * Iz / l^3 * kBendReleaseLef ;
    #   elseif elemReleases(3) == 0 && elemReleases(4) == 1
    #     KbendXY = E * Iz / l^3 * kBendReleaseRig ;
    #   else
    #     KbendXY = zeros(4,4) ;
    #   end

    # 	% bending XZ
    # 	RXYXZ = eye(4) ; RXYXZ(2,2) = -1; RXYXZ(4,4) = -1;
    # 	if     elemReleases(1) == 0 && elemReleases(2) == 0
    # 		KbendXZ = E * Iy / l^3 * RXYXZ * kBendNoRelease * RXYXZ ;
    # 	elseif elemReleases(1) == 1 && elemReleases(2) == 0
    # 		KbendXZ = E * Iy / l^3 * RXYXZ * kBendReleaseLef * RXYXZ ;
    # 	elseif elemReleases(1) == 0 && elemReleases(2) == 1
    # 		KbendXZ = E * Iy / l^3 * RXYXZ * kBendReleaseRig * RXYXZ ;
    # 	else
    # 		KbendXZ = zeros(4,4) ;
    # 	end

    #   Ktorsn = G*J/l * [  1 -1  ; ...
    #                      -1  1  ] ;

    #   KL( LocBendXYdofs , LocBendXYdofs ) = KbendXY ;
    #   KL( LocBendXZdofs , LocBendXZdofs ) = KbendXZ ;
    #   KL( LocTorsndofs  , LocTorsndofs  ) = Ktorsn ;

    #   KGelem = R * KL * R' ;
    #   Finte = KGelem * Ut ;

    #   finteLocalCoor =  R' * Finte ;
    return 1, 1, 1, 1
end

end # module
