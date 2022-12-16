using BenchmarkTools, StaticArrays, LinearAlgebra

# # nodes2dofs

# # actual version
# function nodes2dofs(nodes, degreespernode)
#     n = length(nodes)
#     dofs = zeros(Int, n * degreespernode)
#     for i in (1:n)
#         # print("\ni:", i, "\n")
#         # print((i - 1) * degreespernode .+ Vector(1:degreespernode))
#         # print((degreespernode * (nodes[i] - 1)) .+ Vector(1:degreespernode))
#         dofs[(i-1)*degreespernode.+Vector(1:degreespernode)] = (degreespernode * (nodes[i] - 1)) .+ Vector(1:degreespernode)
#     end
#     return dofs
# end

# # version in materialnonlinearity
# function nodes2dofs_2(nodes, ndofs)
#     n = length(nodes)
#     gdl = zeros(Int64, n * ndofs)
#     vec = Vector(1:ndofs)
#     for i in 1:n
#         gdl[(i-1)*ndofs.+vec] = (nodes[i] - 1) * ndofs .+ vec
#     end
#     return gdl
# end

# # version using StaticArrays
# function nodes2dofs_3(nodes, ndofs)
#     n = length(nodes)
#     gdl = zeros(SizedVector{n * ndofs})
#     vec = SVector{ndofs}(1:ndofs)
#     for i in 1:n
#         gdl[(i-1)*ndofs.+vec] = (nodes[i] - 1) * ndofs .+ vec
#     end
#     return gdl
# end

# function nodes2dofs_4(nodes, ndofs)
#     n = length(nodes)
#     gdl = zeros(MVector{n * ndofs})
#     vec = SVector{ndofs}(1:ndofs)
#     for i in 1:n
#         gdl[(i-1)*ndofs.+vec] = (nodes[i] - 1) * ndofs .+ vec
#     end
#     return gdl
# end

# function nodes2dofs_5(nodes, ndofs)
#     gdl = reduce(vcat, [collect((nodes[i] - 1) * ndofs .+ (1:ndofs)) for i = 1:n])
#     return gdl
# end

# function nodes2dofs_6(nodes, ndofs)
#     gdl = reduce(vcat, [(nodes[i] - 1) * ndofs .+ (1:ndofs) for i = 1:n])
#     return gdl
# end

# nodes = [1, 2, 3, 4]
# ndofs = 6

# @btime nodes2dofs($nodes, $ndofs)
# @btime nodes2dofs_2($nodes, $ndofs)
# @btime nodes2dofs_3($nodes, $ndofs)
# @btime nodes2dofs_4($nodes, $ndofs)
# @btime nodes2dofs_5($nodes, $ndofs)
# @btime nodes2dofs_6($nodes, $ndofs)

# frame geometry

"""
function to linear truss element
"""
function linear_truss(material, geometry, nodalCoords, u)

    # E = material.constitutive_params[1]
    # A = geometry.area
    E = 1.0
    A = 1.0

    diff = nodalCoords[2, :] - nodalCoords[1, :]

    length = sqrt(diff' * diff)
    # println(length)
    c = diff[1] / length
    s = diff[3] / length

    Qloc2glo = [c -s 0 0
        s c 0 0
        0 0 c -s
        0 0 s c]

    Kloc = E * A / length * [1 0 -1 0
        0 0 0 0
        -1 0 1 0
        0 0 0 0]

    Kglo = Qloc2glo * Kloc * transpose(Qloc2glo)

    fint = Kglo * u

    force_vectors = [fint]
    tangent_matrices = [Kglo]

    return force_vectors, tangent_matrices

end

function linear_truss2(material, geometry, nodalCoords::Matrix, u::Vector)

    # E = material.constitutive_params[1]
    # A = geometry.area
    E = 1.0
    A = 1.0

    diff = view(nodalCoords, 2, :) - view(nodalCoords, 1, :)

    length = sqrt(diff' * diff)
    # c = diff[1] / length
    # s = diff[3] / length
    (c, s) = (diff[1], diff[3]) ./ length

    Qloc2glo = [c -s 0 0
        s c 0 0
        0 0 c -s
        0 0 s c]

    Kloc = E * A / length * [1 0 -1 0
        0 0 0 0
        -1 0 1 0
        0 0 0 0]

    Kglo = Qloc2glo * Kloc * Qloc2glo'

    fint = Kglo * u

    force_vectors = [fint]
    tangent_matrices = [Kglo]

    return force_vectors, tangent_matrices

end


function linear_truss3(material, geometry, nodalCoords::Matrix, u::Vector)

    # E = material.constitutive_params[1]
    # A = geometry.area
    E = 1.0
    A = 1.0

    diff = SVector{3}(nodalCoords[2, :]) - SVector{3}(nodalCoords[1, :])

    length = sqrt(diff' * diff)
    # c = diff[1] / length
    # s = diff[3] / length
    (c, s) = (diff[1], diff[3]) ./ length

    Qloc2glo = @SMatrix [c -s 0 0
        s c 0 0
        0 0 c -s
        0 0 s c]

    Kloc = E * A / length * @SMatrix [1 0 -1 0
        0 0 0 0
        -1 0 1 0
        0 0 0 0]

    Kglo = Qloc2glo * Kloc * Qloc2glo'

    fint = Kglo * u

    force_vectors = [fint ]
    tangent_matrices = [Kglo ]

    return force_vectors, tangent_matrices

end

function mats(c, s)
    Q = [c -s 0 0
        s c 0 0
        0 0 c -s
        0 0 s c]
    return Q
end

u = rand(4)
nodalCoords = [[1 2 3]; [3 4 5]]

@btime linear_truss($nothing, $nothing, $nodalCoords, $u)
@btime linear_truss2($nothing, $nothing, $nodalCoords, $u)
@btime linear_truss3($nothing, $nothing, $nodalCoords, $u)

a, b = linear_truss(nothing, nothing, nodalCoords, u)
a3, b3 = linear_truss3(nothing, nothing, nodalCoords, u)
