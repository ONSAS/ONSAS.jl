using ..Materials: SVK
using ..Elements: AbstractElement, AbstractNode, Node, _DEFAULT_LABEL
using ..CrossSections: AbstractCrossSection, area
using ..Utils: eye, row_vector

export Truss

"""
A `Truss` represents a 2D element that transmits axial force only.
### Fields:
- `nodes`          -- stores truss nodes.
- `cross_sections` -- stores the truss cross-section properties.
- `label`          -- stores the truss label.
"""
struct Truss{dim,G<:AbstractCrossSection,T<:Real} <: AbstractElement{dim,T}
    nodes::Vector{<:AbstractNode{dim,T}}
    cross_section::G
    label::Symbol
    function Truss(nodes::Vector{<:AbstractNode{dim,T}}, g::G, label=_DEFAULT_LABEL) where
    {dim,G<:AbstractCrossSection,T<:Real}
        setindex!(nodes[1], 1)
        setindex!(nodes[2], 2)
        new{dim,G,T}(nodes, g, Symbol(label))
    end
end
#TODO: Implement a more general strain function strain(::Truss{E::GreenLagrange})

_strain(l_ini::Number, l_def::Number) = 0.5 * (l_def^2 - l_ini^2) / l_ini^2 # green lagrange strain

function strain(e::Truss, u_e::AbstractVector)
    l_def, l_ini = lengths(e, u_e)
    ϵ = _strain(l_ini, l_def)
end

function internal_forces(m::SVK, e::Truss{dim}, u_e::AbstractVector) where {dim}

    E = m.E
    A = area(cross_section(e))

    u_loc_dofs = view(u_e, 1:2:length(u_e))

    X_ref, X_def = _X_rows(e, u_loc_dofs)
    l_ref, l_def = _lengths(X_ref, X_def, dim)
    _, G = _aux_matrices(dim)
    b_ref, b_def = _aux_b(X_ref, X_def, u_loc_dofs, G, dim)

    ϵ_e = _strain(l_ref, l_def)
    σ_e = E * ϵ_e
    aux = (b_ref + b_def)
    fᵢₙₜ_e = A * σ_e * l_ref * aux'
    Kᵢₙₜ_e = A * σ_e / l_ref * G + E * A * l_ref * (aux' * aux)

    return fᵢₙₜ_e, Kᵢₙₜ_e, σ_e, ϵ_e
end

"Returns local dofs of a truss element."
local_dofs(::Truss{1}) = [Dof(1), Dof(2)]
local_dofs(::Truss{2}) = [Dof(1), Dof(3)]
local_dofs(::Truss{3}) = [Dof(1), Dof(3), Dof(5)]

function _aux_matrices(dim::Integer)
    Bdif = hcat(-eye(dim), eye(dim))
    Ge = Bdif' * Bdif
    return Bdif, Ge
end

function _X_rows(e::Truss{dim}, u_e::AbstractVector) where {dim}
    X_ref_row = coordinates(e)
    X_def_row = X_ref_row + u_e
    return X_ref_row, X_def_row
end

function _aux_b(X_ref_row::AbstractVector, X_def_row::AbstractVector, u_loc_dofs::AbstractVector, G::AbstractMatrix, dim::Integer)

    l_ref, l_def = _lengths(X_ref_row, X_def_row, dim)

    b_ref = 1 / (l_ref^2) * X_ref_row' * G
    b_def = 1 / (l_def^2) * u_loc_dofs' * G

    return b_ref, b_def
end

"Returns deformed and reference lengths of a truss element."
function _lengths(X_ref_row::AbstractVector, X_def_row::AbstractVector, dim::Integer)
    Bdif, _ = _aux_matrices(dim)
    l_ref = sqrt(sum((Bdif * X_ref_row) .^ 2))
    l_def = sqrt(sum((Bdif * X_def_row) .^ 2))

    return l_ref, l_def
end