using ONSAS
using LinearAlgebra

#
# Idea: implementar internal_forces de Frame.jl
#

function ONSAS.internal_forces(m::IsotropicLinearElastic, f::Frame, u_e::AbstractVector)
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
    Izz = CrossSections.Izz(S)
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

# Ejemplo
# ========
m = IsotropicLinearElastic(2.1e11, 0.3)
f = Frame(Node(0.0, 0.0, 0.0), Node(0.75, 0.0, 0.0), Rectangle(0.2, 0.3))
u_e = zeros(12)
# [u1_1, u2_1, u3_1, u1_2, u2_2, u3_2, t1_1, t2_1, t3_1, t1_2, t2_2, t3_2]

internal_forces(m, f, u_e)
