#############################
# Materials interface tests #
#############################
using Test: @testset, @test
using ONSAS.Materials

const TOLERANCE = 1e-2

@testset "ONSAS.Materials SVK" begin

    # random tpyes for each SVK parameter
    types = [Float64, Int, Float32, Int64]
    type1 = rand(types)
    type2 = rand(types)

    E = type1(rand(3e3:10e9)[1])
    ν = type2(rand(0:0.5)[1])

    # SVK for static analysis
    svk_static = SVK(E, ν)
    @test svk_static.E == E
    @test svk_static.ν == ν
    @test svk_static.ρ == nothing

    # SVK for dynamic analysis
    ρ = rand(type1, 1)[1]
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
    @test model(svk_lame) == "SVK"
    @test parameters(svk_lame)[1] ≈ E atol = TOLERANCE
    @test parameters(svk_lame)[2] ≈ ν atol = TOLERANCE
    @test label(svk_lame) == String(label_lame)

end