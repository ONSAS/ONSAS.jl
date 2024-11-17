using Test, LinearAlgebra, ForwardDiff

# Modules to test
using ONSAS.Materials
using ONSAS.LinearElasticMaterials
using ONSAS.IsotropicLinearElasticMaterial
using ONSAS.HyperElasticMaterials
using ONSAS.SVKMaterial
using ONSAS.NeoHookeanMaterial
using ONSAS.HyperElasticMaterial

using ONSAS.Utils

const RTOL = 1e-3

# Steel
E = 210e9
ν = 0.3
G = E / (2 * (1 + ν))
λ = E * ν / ((1 + ν) * (1 - 2 * ν))
K = E / (3 * (1 - 2 * ν))
ρ = 7500.0
mat_label = "steel"

@testset "ONSAS.IsotropicLinearElasticMaterial" begin
    linear_steel_no_density = IsotropicLinearElastic(E, ν)

    @test elasticity_modulus(linear_steel_no_density) == E
    @test poisson_ratio(linear_steel_no_density) == ν
    @test shear_modulus(linear_steel_no_density) == G
    @test bulk_modulus(linear_steel_no_density) == K
    @test isnothing(density(linear_steel_no_density))
    @test label(linear_steel_no_density) == NO_LABEL

    linear_steel = IsotropicLinearElastic(; λ = λ, G = G, ρ = ρ, label = mat_label)
    @test label(linear_steel) == Symbol(mat_label)
    @test density(linear_steel) == ρ
    @test elasticity_modulus(linear_steel)≈E rtol=RTOL
    @test poisson_ratio(linear_steel)≈ν rtol=RTOL
    @test shear_modulus(linear_steel)≈G rtol=RTOL
    @test bulk_modulus(linear_steel)≈K rtol=RTOL

    # Test constitutive driver
    ϵᵢ = 0.18375
    ϵⱼ = 0.435
    ϵᵏ = 1.14
    γᵢⱼ = 0.2925
    γⱼₖ = 0.72
    γₖᵢ = 0.51

    # Consitutive tensor
    𝐶 = [λ+2G λ λ 0 0 0
         λ λ+2G λ 0 0 0
         λ λ λ+2G 0 0 0
         0 0 0 G 0 0
         0 0 0 0 G 0
         0 0 0 0 0 G]

    ϵ = Symmetric([ϵᵢ γᵢⱼ γₖᵢ
                   γᵢⱼ ϵⱼ γⱼₖ
                   γₖᵢ γⱼₖ ϵᵏ])

    ϵ_vec = voigt(ϵ, 2)
    σ_vogit = 𝐶 * ϵ_vec
    σ_expected = Symmetric([σ_vogit[1] σ_vogit[6] σ_vogit[5]
                            σ_vogit[6] σ_vogit[2] σ_vogit[4]
                            σ_vogit[5] σ_vogit[4] σ_vogit[3]])

    σ = Symmetric(zeros(3, 3))
    ∂σ∂ϵ = zeros(6, 6)

    σ, ∂σ∂ϵ = stress!(σ, ∂σ∂ϵ, linear_steel, ϵ)

    @test σ≈σ_expected rtol=RTOL
end

# More soft hyperelastic material
Ghyper = μ = 0.3846
λhyper = 0.5769
Khyper = λhyper + 2 * Ghyper / 3
# Green-Lagrange strain tensor for testing
𝔼 = Symmetric([0.18375 0.2925 0.51
               0.2925 0.435 0.72
               0.51 0.72 1.14])

@testset "ONSAS.SVKMaterial + ONSAS.HyperElasticMaterial" begin

    # SVK for static analysis
    svk_static = SVK(λ, G)

    @test lame_parameters(svk_static) == (λ, G)
    @test isnothing(density(svk_static))
    @test elasticity_modulus(svk_static)≈E rtol=RTOL
    @test shear_modulus(svk_static) == G
    @test bulk_modulus(svk_static)≈K rtol=RTOL
    @test poisson_ratio(svk_static) == ν

    # SVK for dynamic analysis
    svk_dynamic = SVK(; E = E, ν = ν, ρ = ρ, label = mat_label)
    @test density(svk_dynamic) == ρ
    @test collect(lame_parameters(svk_dynamic))≈[λ, G] rtol=RTOL
    @test label(svk_dynamic) == Symbol(mat_label)
    # SVK strain energy
    strain_energy_svk(𝔼, λ::Real, G::Real) = (λ / 2) * tr(𝔼)^2 + G * tr(𝔼^2)
    @test strain_energy(svk_dynamic, 𝔼) ==
          strain_energy_svk(𝔼, lame_parameters(svk_static)...)

    l = "svk_HyperElastic"
    svk_hyper = HyperElastic([λhyper, Ghyper], strain_energy_svk, l)
    svk = SVK(λhyper, Ghyper)

    𝕊_test = Symmetric([1.15596 0.224991 0.392292
                        0.224991 1.34922 0.553824
                        0.392292 0.553824 1.89151])

    ∂𝕊∂𝔼_test = [1.3461 0.5769 0.5769 0.0 0.0 0.0
                 0.5769 1.3461 0.5769 0.0 0.0 0.0
                 0.5769 0.5769 1.3461 0.0 0.0 0.0
                 0.0 0.0 0.0 0.3846 0.0 0.0
                 0.0 0.0 0.0 0.0 0.3846 0.0
                 0.0 0.0 0.0 0.0 0.0 0.3846]

    @test parameters(svk_hyper) == [λhyper, Ghyper]
    @test isnothing(density(svk_hyper))
    @test label(svk_hyper) == Symbol(l)

    # Constitutive driver svk type SVK
    𝕊_svk = Symmetric(zeros(3, 3))
    ∂𝕊∂𝔼_svk = zeros(6, 6)
    cosserat_stress!(𝕊_svk, ∂𝕊∂𝔼_svk, svk, 𝔼)
    @test 𝕊_svk≈𝕊_test rtol=RTOL
    @test ∂𝕊∂𝔼_svk≈∂𝕊∂𝔼_test rtol=RTOL
    # Constitutive driver HyperElasticMateiral
    𝕊_hyper = Symmetric(zeros(3, 3))
    ∂𝕊∂𝔼_hyper = zeros(6, 6)
    cosserat_stress!(𝕊_hyper, ∂𝕊∂𝔼_hyper, svk_hyper, 𝔼)
    @test 𝕊_hyper≈𝕊_test rtol=RTOL
    @test ∂𝕊∂𝔼_svk≈∂𝕊∂𝔼_test rtol=RTOL
end

@testset "ONSAS.SVKMaterial + ONSAS.NeoHookeanMaterial" begin
    neo = NeoHookean(K, G)
    @test bulk_modulus(neo) == K
    @test shear_modulus(neo) == G
    @test collect(lame_parameters(neo))≈[λ, G] rtol=RTOL
    @test poisson_ratio(neo)≈ν rtol=RTOL
    @test elasticity_modulus(neo)≈E rtol=RTOL

    # NeoHookean defined with ρ E and  ν
    neo_withρ = NeoHookean(; E = E, ν = ν, ρ = ρ, label = mat_label)
    @test bulk_modulus(neo_withρ)≈K rtol=RTOL
    @test shear_modulus(neo_withρ)≈G rtol=RTOL
    @test collect(lame_parameters(neo_withρ))≈[λ, G] rtol=RTOL
    @test poisson_ratio(neo_withρ)≈ν rtol=RTOL
    @test elasticity_modulus(neo_withρ)≈E rtol=RTOL
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
        return Ψ = μ / 2 * (I₁ - 2 * log(J)) + K / 2 * (J - 1)^2
    end

    neo_hyper = HyperElastic([bulk_modulus(neo_flexible), shear_modulus(neo_flexible)],
        strain_energy_neo, l)

    𝕊_hyper = Symmetric(zeros(3, 3))
    ∂𝕊∂𝔼_hyper = zeros(6, 6)
    𝕊_neo = Symmetric(zeros(3, 3))
    ∂𝕊∂𝔼_neo = zeros(6, 6)

    cosserat_stress!(𝕊_hyper, ∂𝕊∂𝔼_hyper, neo_hyper, 𝔼)
    cosserat_stress!(𝕊_neo, ∂𝕊∂𝔼_neo, neo_flexible, 𝔼)

    @test 𝕊_hyper≈𝕊_neo rtol=RTOL
    @test ∂𝕊∂𝔼_hyper≈∂𝕊∂𝔼_neo rtol=RTOL
end
