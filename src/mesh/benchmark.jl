using BenchmarkTools, StaticArrays, LinearAlgebra

# actual version
function nodes2dofs(nodes, degreespernode)
    n = length(nodes)
    dofs = zeros(Int, n * degreespernode)
    for i in (1:n)
        # print("\ni:", i, "\n")
        # print((i - 1) * degreespernode .+ Vector(1:degreespernode))
        # print((degreespernode * (nodes[i] - 1)) .+ Vector(1:degreespernode))
        dofs[(i-1)*degreespernode.+Vector(1:degreespernode)] = (degreespernode * (nodes[i] - 1)) .+ Vector(1:degreespernode)
    end
    return dofs
end

# version in materialnonlinearity
function nodes2dofs_2(nodes, ndofs)
    n = length(nodes)
    gdl = zeros(Int64, n * ndofs)
    vec = Vector(1:ndofs)
    for i in 1:n
        gdl[(i-1)*ndofs.+vec] = (nodes[i] - 1) * ndofs .+ vec
    end
    return gdl
end

# version using StaticArrays
function nodes2dofs_3(nodes, ndofs)
    n = length(nodes)
    gdl = zeros(SizedVector{n * ndofs})
    vec = SVector{ndofs}(1:ndofs)
    for i in 1:n
        gdl[(i-1)*ndofs.+vec] = (nodes[i] - 1) * ndofs .+ vec
    end
    return gdl
end

function nodes2dofs_4(nodes, ndofs)
    n = length(nodes)
    gdl = zeros(MVector{n * ndofs})
    vec = SVector{ndofs}(1:ndofs)
    for i in 1:n
        gdl[(i-1)*ndofs.+vec] = (nodes[i] - 1) * ndofs .+ vec
    end
    return gdl
end

function nodes2dofs_5(nodes, ndofs)
    n = length(nodes)
    gdl = reduce(vcat, [collect((nodes[i] - 1) * ndofs .+ (1:ndofs)) for i = 1:n])
    return gdl
end

function nodes2dofs_6(nodes, ndofs)
    n = length(nodes)
    gdl = reduce(vcat, [(nodes[i] - 1) * ndofs .+ (1:ndofs) for i = 1:n])
    return gdl
end

function nodes2dofs_7(nodes::Vector{Int}, ndofs::Int)
    n = length(nodes)
    gdl = Vector{Int}(undef, n * ndofs)
    @inbounds for i in 1:n
        α = (nodes[i] - 1) * ndofs
        for j in 1:ndofs
            β = (i - 1) * ndofs
            gdl[β+j] = α + j
        end
    end
    return gdl
end

nodes = collect((1:300))
ndofs = 6

@btime nodes2dofs($nodes, $ndofs)
@btime nodes2dofs_2($nodes, $ndofs)
@btime nodes2dofs_3($nodes, $ndofs)
@btime nodes2dofs_4($nodes, $ndofs)
@btime nodes2dofs_5($nodes, $ndofs)
@btime nodes2dofs_6($nodes, $ndofs)
@btime nodes2dofs_7($nodes, $ndofs)

