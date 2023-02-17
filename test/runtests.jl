using Test
using SafeTestsets: @safetestset

MODULES = ["materials.jl", "cross_sections.jl"]
EXAMPLES = [
# joinpath("..", "examples", "vonMisesTruss", "von_misses_truss.jl"),
#joinpath("..", "examples", "uniaxialExtension", "uniaxialExtension.jl")
]

function test(files::Vector)
    @testset "ONSAS.jl" begin
        foreach(test, files)
    end
end

function test(file::String)
    @info "Testing $file..."
    path = joinpath(@__DIR__, file)
    @eval @time @safetestset $file begin
        include($path)
    end
end

test([MODULES; EXAMPLES]);
