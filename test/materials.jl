#############################
# Materials interface tests #
#############################
using Test: @testset, @test
using ONSAS.Materials

const TOLERANCE = 1e-2

@testset "ONSAS.Materials SVK" begin

    E = 2e9
    ν = 1 / 3

    # SVK for static analysis
    svk_static = SVK(E, ν)
    @test svk_static.E == E
    @test svk_static.ν == ν
    @test svk_static.ρ == nothing

    # SVK for dynamic analysis
    ρ = 1454
    svk_dynamic = SVK(E, ν, ρ, "mat1")
    @test svk_dynamic.E == E
    @test svk_dynamic.ν == ν
    @test svk_dynamic.ρ == ρ

    # SVK with lamé parameters (λ, G)
    G = E / (2 * (1 + ν))
    λ = E * ν / ((1 + ν) * (1 - 2ν))
    label_lame = :test_label
    svk_lame = SVK(λ=λ, G=G, ρ=ρ, label=label_lame)

    # Compute lamé parameters back
    lame_params = lame_parameters(svk_lame)
    @test lame_params[1] ≈ λ atol = TOLERANCE
    @test lame_params[2] ≈ G atol = TOLERANCE

    # Test Abstract Material interface
    @test parameters(svk_lame)[1] ≈ E atol = TOLERANCE
    @test parameters(svk_lame)[2] ≈ ν atol = TOLERANCE
    @test label(svk_lame) == label_lame

end