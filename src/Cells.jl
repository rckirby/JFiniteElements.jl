module Cells

export Cell, Simplex, UFCSimplex, val2val, getSpatialDimension, getVertexCoords

abstract Cell

abstract Simplex{T<:Val} <: Cell

immutable UFCSimplex{T<:Val} <: Simplex{T} end

immutable TensorProductCell{C1<:Cell, C2<:Cell} <: Cell end


# helper functions used in UFC implementation
function comp(d, vals)
    x = filter((i)->!(i in vals), 1:d)
    return ntuple((i)->x[i], d-length(vals))
end

function concat(a::Tuple, b::Tuple)
    clist = vcat([ai for ai in a],
                 [bi for bi in b])
    return ntuple((i)->clist[i], length(clist),)
end

val2val{d}(::Type{Val{d}}) = d

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

# This is how to get spatial dimension to work on subtypes, defining it
# only on the leaf node in the hierarchy where we branch
# I have a "rootleaf" that will say each kind of simplex is in fact a simplex,
# then just provide the spatial dimension on "simplex".  I can do the same for
# other kinds, I hope.
getSpatialDimension{T}(::Type{Simplex{T}}) = val2val(T)
getSpatialDimension{S}(::Type{S}) = getSpatialDimension(rootleaf(S))
function getSpatialDimension{C1,C2}(::Type{TensorProductCell{C1,C2}})
    return getSpatialDimension(C1) + getSpatialDimension(C2)
end

rootleaf{T}(::Type{UFCSimplex{T}}) = Simplex{T}

getNumVertices{T}(::Type{Simplex{T}}) = 1 + val2val(T)
getNumVertices{S}(::Type{S}) = getNumVertices(rootleaf(S))
function getNumVertices{C1,C2}(::TensorProductCell{C1,C2})
    return getNumVertices(C1) * getNumVertices(C2)
end


# this is generated since it's once per type, so can
# compute array at compile-time and return it at run-time.
@generated function getVertexCoords{T}(::Type{UFCSimplex{T}})
    d = val2val(T)
    x = zeros(d, d+1)
    for i=1:d
        x[i,i+1] = 1.0
    end
    return :($x)
end

@generated function getCellTopology{T}(::Type{UFCSimplex{T}})
    d = val2val(T)
    
    D = Dict()
    D[0] = [(i,) for i=1:d+1]
    
    # Q: Can we meta-program this to get the right number of
    # of Tuple{Int,...} as the type of D[d-dim]?
    for dim=1:(d-1)
        off_verts = get_subsequences(1, d+1, dim)
        D[d-dim] = [comp(d+1, ov) for ov in off_verts]
    end

    # cell itself is pretty easy
    D[d] = [ntuple((i)->i, d+1)]
    return :($D)
end

@generated function getVertexCoords{C1,C2}(::Type{TensorProductCell{C1,C2}})
    coords1 = getVertexCoords(C1)
    coords2 = getVertexCoords(C2)

    sdim1 = getSpatialDimension(C1)
    sdim2 = getSpatialDimension(C2)
    sdim = sdim1 + sdim2
    verts1 = getVertexCoords(C1)
    verts2 = getVertexCoords(C2)
    nverts1 = getNumVertices(C1)
    nverts2 = getNumVertices(C2)

    nverts = nverts1 * nverts2

    verts = zeros(sdim, nverts)
    vert_cur = 1
    for i=1:nverts2
        verts[(sdim1+1):(sdim1+sdim2),vert_cur] = verts2[:,i]
        for j=1:nverts1
            verts[1:sdim1,vert_cur] = verts1[:,j]
            vert_cur += 1
        end
    end

    return :($verts)
end

end
