##########################
# Materials module tests #
##########################
using Test: @testset, @test
using ONSAS.Materials
using ONSAS.Utils: eye
using LinearAlgebra: Symmetric, tr, det, inv
const RTOL = 1e-3

# Steel 
E = 210e9
ν = 0.3
G = E / (2 * (1 + ν))
λ = E * ν / ((1 + ν) * (1 - 2 * ν))
K = E / (3 * (1 - 2 * ν))
ρ = 7500.0
mat_label = "steel"

# More soft hyperelastic material   
Ghyper = μ = 0.3846
λhyper = 0.5769
Khyper = λhyper + 2 * Ghyper / 3
# Green-Lagrange strain tensor for testing
𝔼 = Symmetric(
    [
        0.18375 0.2925 0.51
        0.2925 0.435 0.72
        0.51 0.72 1.14
    ]
)


@testset "ONSAS.Materials.SVK" begin

    # SVK for static analysis
    svk_static = SVK(λ, G)

    @test lame_parameters(svk_static) == (λ, G)
    @test density(svk_static) == nothing
    @test elasticity_modulus(svk_static) ≈ E rtol = RTOL
    @test shear_modulus(svk_static) == G
    @test bulk_modulus(svk_static) ≈ K rtol = RTOL
    @test poisson_ratio(svk_static) == ν

    # SVK for dynamic analysis
    svk_dynamic = SVK(E=E, ν=ν, ρ=ρ, label=mat_label)
    @test density(svk_dynamic) == ρ
    @test lame_parameters(svk_dynamic) |> collect ≈ [λ, G] rtol = RTOL
    @test label(svk_dynamic) == Symbol(mat_label)
    # SVK strain energy
    strain_energy_svk(𝔼, λ::Real, G::Real) = (λ / 2) * tr(𝔼)^2 + G * tr(𝔼^2)
    @test strain_energy(svk_dynamic, 𝔼) == strain_energy_svk(𝔼, lame_parameters(svk_static)...)


    l = "svk_HyperElastic"
    svk_hyper = HyperElastic([λhyper, Ghyper], strain_energy_svk, l)
    svk = SVK(λhyper, Ghyper)


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


    @test parameters(svk_hyper) == [λhyper, Ghyper]
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

@testset "ONSAS.Materials.NeoHookean" begin

    neo = NeoHookean(K, G)
    @test bulk_modulus(neo) == K
    @test shear_modulus(neo) == G
    @test lame_parameters(neo) |> collect ≈ [λ, G] rtol = RTOL
    @test poisson_ratio(neo) ≈ ν rtol = RTOL
    @test elasticity_modulus(neo) ≈ E rtol = RTOL

    # NeoHookean defined with ρ E and  ν
    neo_withρ = NeoHookean(E=E, ν=ν, ρ=ρ, label=mat_label)
    @test bulk_modulus(neo_withρ) ≈ K rtol = RTOL
    @test shear_modulus(neo_withρ) ≈ G rtol = RTOL
    @test lame_parameters(neo_withρ) |> collect ≈ [λ, G] rtol = RTOL
    @test poisson_ratio(neo_withρ) ≈ ν rtol = RTOL
    @test elasticity_modulus(neo_withρ) ≈ E rtol = RTOL
    @test label(neo_withρ) == Symbol(mat_label)

    # More flexible noe-hookean to test strain and stresses
    neo_flexible = NeoHookean(Khyper, Ghyper)

    # Create an hyper-elastic material with the same strain energy and test 𝕊 and 𝔼
    l = "neo_HyperElastic"
    function strain_energy_neo(𝔼::AbstractMatrix, K::Real, μ::Real)
        # Right hand Cauchy strain tensor
        ℂ = Symmetric(2 * 𝔼 + eye(3))
        J = sqrt(det(ℂ))
        # First invariant
        I₁ = tr(ℂ)
        # Strain energy function 
        Ψ = μ / 2 * (I₁ - 2 * log(J)) + K / 2 * (J - 1)^2
    end

    neo_hyper = HyperElastic(
        [bulk_modulus(neo_flexible), shear_modulus(neo_flexible)], strain_energy_neo, l
    )

    𝕊_hyper, ∂𝕊∂𝔼_hyper = cosserat(neo_hyper, 𝔼)
    𝕊_neo, ∂𝕊∂𝔼_neo = cosserat(neo_flexible, 𝔼)

    @test 𝕊_hyper ≈ 𝕊_neo rtol = RTOL
    @test ∂𝕊∂𝔼_hyper ≈ ∂𝕊∂𝔼_neo rtol = RTOL

end