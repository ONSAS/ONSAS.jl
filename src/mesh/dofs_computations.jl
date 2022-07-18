"""
computeElemsByMEBI
This function computes the vector of the element indexes for a given vector of MEBI parameter
"""
function computeElemsByMEBI( MEBIVals, MEBIVec )

    # get the relevant BCs definitions
    MEBIUnique = sort( unique( MEBIVals ) )     # get and sort them
    MEBIUnique[1]==0 && popfirst!( MEBIUnique ) # and remove zeros

    elements = [] 
    # loop
    print(MEBIUnique)
    for MEBIval in MEBIUnique
        print("\nMEBI:", MEBIval, "\n" )
        # find the elements with the current MEBIval
        aux = findall( x->x==MEBIval, MEBIVec )
        
        # if there are any, add them to the list
        length(aux) > 0 && push!( elements, aux )
    end
    return elements
end


"""
nodes2dofs computes the vector of dofs for an input vector of nodes
"""
function nodes2dofs( nodes , degreespernode )
    n    = length(nodes);
    dofs = zeros( Int, n*degreespernode ) ;
    for i in (1:n);
        print("\ni:",i,"\n")
        print( (i-1)*degreespernode .+ Vector(1:degreespernode) )
        print((degreespernode*(nodes[i]-1)).+Vector(1:degreespernode))
        dofs[ (i-1)*degreespernode .+ Vector(1:degreespernode) ] = (degreespernode*(nodes[i]-1)).+Vector(1:degreespernode)  ;
    end
    return dofs;
end



function compute_sparse_indexes( matrix, dofsElem )

    num_cols = size( matrix, 2)
    num_dofs = length( dofsElem )

    num_cols != num_cols && error("number of dofs must be equal to matrix cols")

    row_indexes = Vector{Int}()
    col_indexes = Vector{Int}()
    vals_vector = Vector{Float64}()

    for i in (1:num_cols)
        
        append!(row_indexes, dofsElem )
        append!(col_indexes, dofsElem[i]*ones(num_dofs) )
        append!(vals_vector, matrix[:,i])
    end

    return row_indexes, col_indexes, vals_vector
end













