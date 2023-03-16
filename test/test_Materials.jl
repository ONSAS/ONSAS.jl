##########################
# Materials module tests #
##########################
using Test: @testset, @test
using ONSAS.Materials
using LinearAlgebra: Symmetric, tr
const RTOL = 1e-3

strain_energy_svk(𝔼, λ::Real, G::Real) = (λ / 2) * tr(𝔼)^2 + G * tr(𝔼^2)

@testset "ONSAS.Materials.SVK" begin

    # Steel 
    E = 210e9
    ν = 0.3
    G = E / (2 * (1 + ν))
    λ = E * ν / ((1 + ν) * (1 - 2 * ν))
    K = E / (3 * (1 - 2 * ν))

    # SVK for static analysis
    svk_static = SVK(λ, G)

    @test lame_parameters(svk_static) == (λ, G)
    @test density(svk_static) == nothing
    @test elasticity_modulus(svk_static) ≈ E rtol = RTOL
    @test shear_modulus(svk_static) == G
    @test bulk_modulus(svk_static) ≈ K rtol = RTOL
    @test poisson_ratio(svk_static) == ν

    # SVK for dynamic analysis
    ρ = 7500.0
    label_lame = "steel"
    svk_dynamic = SVK(E=E, ν=ν, ρ=ρ, label=label_lame)
    @test density(svk_dynamic) == ρ
    @test lame_parameters(svk_dynamic)[1] ≈ λ rtol = RTOL
    @test lame_parameters(svk_dynamic)[2] ≈ G rtol = RTOL

    𝔼 = Symmetric(rand(3, 3))

    @test strain_energy(svk_dynamic, 𝔼) == strain_energy_svk(𝔼, lame_parameters(svk_static)...)
    @test label(svk_dynamic) == Symbol(label_lame)

end


@testset "ONSAS.Materials.HyperElastic(SVK)" begin

    λ = 0.5769
    G = μ = 0.3846

    𝕊_test = Symmetric(
        [
            1.15596 0.224991 0.392292
            0.224991 1.34922 0.553824
            0.392292 0.553824 1.89151
        ]
    )

    ∂𝕊∂𝔼_test = [
        1.3461 0.5769 0.5769 0.0 0.0 0.0
        0.5769 1.3461 0.5769 0.0 0.0 0.0
        0.5769 0.5769 1.3461 0.0 0.0 0.0
        0.0 0.0 0.0 0.3846 0.0 0.0
        0.0 0.0 0.0 0.0 0.3846 0.0
        0.0 0.0 0.0 0.0 0.0 0.3846
    ]


    𝔼 = Symmetric(
        [
            0.18375 0.2925 0.51
            0.2925 0.435 0.72
            0.51 0.72 1.14
        ]
    )

    svk = SVK(λ, G)
    # Create a generic HyperElastic material with an SVK   
    l = "svk_HyperElastic"
    svk_hyper = HyperElastic([λ, G], strain_energy_svk, l)

    @test parameters(svk_hyper) == [λ, G]
    @test density(svk_hyper) == nothing
    @test label(svk_hyper) == Symbol(l)

    # Constitutive driver svk type SVK
    𝕊_svk, ∂𝕊∂𝔼_svk = cosserat(svk, 𝔼)

    @test 𝕊_svk ≈ 𝕊_test rtol = RTOL
    @test ∂𝕊∂𝔼_svk ≈ ∂𝕊∂𝔼_test rtol = RTOL

    𝕊_hyper, ∂𝕊∂𝔼_hyper = cosserat(svk_hyper, 𝔼)

    @test 𝕊_hyper ≈ 𝕊_test rtol = RTOL
    @test ∂𝕊∂𝔼_svk ≈ ∂𝕊∂𝔼_test rtol = RTOL

end