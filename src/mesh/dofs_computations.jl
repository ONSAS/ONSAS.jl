"""
computeElemsByMEBI
This function computes the vector of the element indexes for a given vector of MEBI parameter
"""
function computeElemsByMEBI(MEBIVals, MEBIVec)

    # get the relevant BCs definitions
    MEBIUnique = sort(unique(MEBIVals))     # get and sort them
    MEBIUnique[1] == 0 && popfirst!(MEBIUnique) # and remove zeros

    elements = []
    # loop
    print(MEBIUnique)
    for MEBIval in MEBIUnique
        print("\nMEBI:", MEBIval, "\n")
        # find the elements with the current MEBIval
        aux = findall(x -> x == MEBIval, MEBIVec)

        # if there are any, add them to the list
        length(aux) > 0 && push!(elements, aux)
    end
    return elements
end


"""
nodes2dofs computes the vector of dofs for an input vector of nodes
"""
function nodes2dofs(nodes::Vector{Integer}, ndofs::Int)
    n = length(nodes)
    gdl = Vector{Int}(undef, n * ndofs)
    @inbounds for i in 1:n
        α = (nodes[i] - 1) * ndofs
        β = (i - 1) * ndofs
        for j in 1:ndofs
            gdl[β+j] = α + j
        end
    end
    return gdl
end


function nodes2dofs(nodes::Number, ndofs::Int)
    n = length(nodes)
    gdl = Vector{Int}(undef, n * ndofs)
    @inbounds for i in 1:n
        α = (nodes[i] - 1) * ndofs
        β = (i - 1) * ndofs
        for j in 1:ndofs
            gdl[β+j] = α + j
        end
    end
    return gdl
end





function compute_sparse_indexes(matrix, dofsElem)

    num_cols = size(matrix, 2)
    num_dofs = length(dofsElem)

    num_cols != num_cols && error("number of dofs must be equal to matrix cols")

    row_indexes = Vector{Int}()
    col_indexes = Vector{Int}()
    vals_vector = Vector{Float64}()

    for i in (1:num_cols)

        append!(row_indexes, dofsElem)
        append!(col_indexes, dofsElem[i] * ones(num_dofs))
        append!(vals_vector, matrix[:, i])
    end

    return row_indexes, col_indexes, vals_vector
end













