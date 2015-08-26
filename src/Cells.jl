module Cells

export Cell, Simplex, UFCSimplex

abstract Cell

abstract Simplex <: Cell

# returns complement of vals considered as a subsequence of (1...d)
function comp(d, vals)
    x = filter((i)->!(i in vals), 1:d)
    return ntuple(d-length(vals), (i)->x[i])
end

function concat(a::Tuple, b::Tuple)
    clist = vcat([ai for ai in a],
                 [bi for bi in b])
    return ntuple(length(clist), (i)->clist[i])
end
             

# produces subsequences of a:b of length l,
# sorted in lexicographic order
function get_subsequences(a, b, l)
    if b < a
        return []
    elseif l == 1
        return [(i,) for i=a:b]
    else
        result = Array(Tuple,0)
        for c=a:b
            for x in get_subsequences(c+1, b, l-1)
                push!(result, concat((c,), x))
            end
        end
        return result
    end
end

immutable UFCSimplex{d} <: Simplex
    x::Array{Float64, 2}   # Vertices
    D::Dict                # Topology
    function UFCSimplex()
        x = zeros(d, d+1)
        for i=1:d
            x[i,i+1] = 1.0
        end

        D = Dict()
        D[0] = [(i,) for i=1:d+1]

        # Q: Can we meta-program this to get the right number of
        # of Tuple{Int,...} args as the type of D[d-dim]?
        for dim=1:(d-1)
            off_verts = get_subsequences(1, d+1, dim)
            D[d-dim] = [comp(d+1, ov) for ov in off_verts]
        end

        # cell itself is pretty easy
        D[d] = [ntuple(d+1, (i)->i)]
        new(x, D)
    end
end

end
