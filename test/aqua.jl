using Aqua
using ONSAS

@testset "Aqua.jl" begin
    Aqua.test_all(ONSAS;
        ambiguities = false)
    Aqua.test_ambiguities(ONSAS; broken = true)
end
