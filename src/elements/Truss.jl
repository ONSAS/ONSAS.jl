using ..Materials: SVK
using ..Utils: ScalarWrapper, eye, row_vector

export Truss

"""
A `Truss` represents a 2D element that transmits axial force only.
### Fields:
- `nodes` -- stores truss nodes.
- `material`  -- stores truss material.
- `geometry` -- stores the truss cross-section properties.
"""
struct Truss{dim,M,G} <: AbstractElement{dim,M}
    nodes::Vector{<:AbstractNode{dim}}
    material::M
    geometry::G
    label::ScalarWrapper{Symbol}
end


function Truss(g::G, dim::Integer=3, label=_DEFAULT_LABEL) where {G<:AbstractCrossSection}
    Truss(Vector{Node{dim,Float64}}(), nothing, g, ScalarWrapper(label))
end

function Truss(m::M, g::G, dim::Integer=3, label=_DEFAULT_LABEL) where {M<:AbstractMaterial,G<:AbstractCrossSection}
    Truss(Vector{Node{dim,Float64}}(), m, g, ScalarWrapper(label))
end

element_type(::Truss) = Truss
num_nodes(::Truss) = 2
dofs_per_node(::Truss{1}) = [Dof(:uₓ, 1)]
dofs_per_node(::Truss{2}) = [Dof(:uₓ, 1), Dof(:uⱼ, 3)]
dofs_per_node(::Truss{3}) = [Dof(:uₓ, 1), Dof(:uⱼ, 3), Dof(:uₖ, 5)]

set_material(t::Truss, m::M) where {M<:AbstractMaterial} = element_type(t)(m, geometry(t), dimension(t), label(t))

function _aux_matrices(dim::Integer)
    Bdif = hcat(-eye(dim), eye(dim))
    Ge = Bdif' * Bdif
    return Bdif, Ge
end

function _X_rows(e::Truss{dim}, u_e::AbstractVector) where {dim}
    X_ref_row = row_vector(coordinates(e))
    X_def_row = row_vector(coordinates(e) + u_e)
    return X_ref_row, X_def_row
end

function _aux_b(X_ref_row::AbstractVector, X_def_row::AbstractVector, G::AbstractMatrix, dim::Integer)

    l_ref, l_def = _lengths(X_ref_row, X_def_row, dim)

    b_ref = 1 / (l_ref^2) * X_ref_row' * G
    b_def = 1 / (l_def^2) * X_def_row' * G

    return b_ref, b_def
end

"Returns deformed and reference lengths of a truss element."
function _lengths(X_ref_row::AbstractVector, X_def_row::AbstractVector, dim::Integer)
    Bdif, _ = _aux_matrices(dim)
    l_ref = sqrt(sum((Bdif * X_ref_row) .^ 2))
    l_def = sqrt(sum((Bdif * X_def_row) .^ 2))

    return l_ref, l_def
end

#TODO: Implement a more general strain function strain(::Truss{E::GreenLagrange})
_strain(l_ini::Number, l_def::Number) = 0.5 * (l_def^2 - l_ini^2) / l_ini^2 # green lagrange strain

function strain(e::Truss{dim}, u_e::AbstractVector) where {dim}
    l_def, l_ini = lengths(e, u_e)
    ϵ = _strain(l_ini, l_def)
end

_stress(l_ini::Number, l_def::Number, E::Number) = E * _strain(l_ini, l_def) # green lagrange stress 

function stress(e::Truss{dim,SVK}, u_e::AbstractVector) where {dim}
    l_def, l_ini = lengths(e, u_e)
    σ = _stress(l_ini, l_def, material(e).E)
end


function stiffness_matrix(e::Truss{dim,SVK}, u_e::AbstractVector) where {dim}

    E = material(e).E
    A = area(geometry(e))

    X_ref, X_def = _X_rows(e, u_e)
    l_ref, l_def = _lengths(X_ref, X_def, dim)
    _, G = _aux_matrices(dim)
    b_ref, b_def = _aux_b(X_ref, X_def, G, dim)

    σ = _stress(l_ref, l_def, E)
    Kₑ = A * σ / l_ref * G + E * A * l_ref * ((b_ref + b_def)' * (b_ref + b_def))
end


function internal_force(e::Truss{dim,SVK}, u_e::AbstractVector) where {dim}

    E = material(e).E
    A = area(geometry(e))

    X_ref, X_def = _X_rows(e, u_e)
    l_ref, l_def = _lengths(X_ref, X_def, dim)
    _, G = _aux_matrices(dim)
    b_ref, b_def = _aux_b(X_ref, X_def, G, dim)

    σ = _stress(l_ref, l_def, E)
    fₑ = A * σ * l_ref * (b_ref + b_def)'
end

