##########################
# Materials module tests #
##########################
using Test: @testset, @test
using ONSAS.Materials
const RTOL = 1e-3

@testset "ONSAS.Materials SVK" begin

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

    @test strain_energy(svk_dynamic) == :((λ / 2) * tr(𝔼)^2 + G * tr(𝔼^2))
    @test label(svk_dynamic) == Symbol(label_lame)

end