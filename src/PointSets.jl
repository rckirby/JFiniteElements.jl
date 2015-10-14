module PointSets

using Cells

import Base.*

export EquispacedLattice, Lattice, PointSet, size, points
export TensorProductPointSet, *, BoundaryMidpointSet

abstract PointSet

abstract Lattice{C<:Cell, S<:Val} <: PointSet

immutable EquispacedLattice{C<:Cell, S<:Val} <: Lattice end

immutable CellMidpointSet{C<:Cell} <: PointSet end
immutable BoundaryMidpointSet{C<:Cell} <: PointSet end

immutable TensorProductPointSet{P1<:PointSet, P2<:PointSet} end

rootleaf{S,T}(::Type{EquispacedLattice{S,T}}) = Lattice{S,T}
size{T}(::Type{T}) = size(rootleaf(T))

type2type{T}(::Type{Type{T}}) = T

cellType{C<:Cell,S<:Val}(::Type{Lattice{C,S}}) = C
cellType{T}(::Type{T}) = cellType(rootleaf(T))


size{C<:Cell}(::Type{CellMidpointSet{C}}) = 1
cellType{C<:Cell}(::Type{CellMidpointSet{C}}) = C
@generated function points{C<:Cell}(::Type{CellMidpointSet{C}})
    verts = getVertexCoords(C)
    # average value of vertices is the midpoint.  Always true?
    v0 = reshape(sum(verts, 2) / size(verts,2))
    return :($v0)
end

@generated function size{C<:Cell}(::Type{BoundaryMidpointSet{C}})
    D = getCellTopology(C)
    sd = spatialDimension(C)
    sz = len(D[sd-1])
    return :($sz)
end

@generated function points{C<:Cell}(::Type{BoundaryMidpointSet{C}})
    verts = getVertexCoords(C)
    D = getCellTopology(C)
    sd = spatialDimension(C)
    num_bd_facets = length(D[sd-1])
    bmps = zeros(sd, num_bd_facets)
    for k=1:num_bd_facets
        vert_ids = D[sd-1][k]
        facet_size = length(vert_ids)
        for j=1:facet_size
            bmps[:,k] += verts[:, vert_ids[j]]
        end
        bmps[:,k] /= facet_size
    end
    return :($bmps)
    
end

@generated function size{S<:Simplex,T<:Val}(::Type{Lattice{S,T}})
    dim = spatialDimension(S)
    n = val2val(T)
    np = 1
    for d=1:dim
        np *= (n+d)
    end
    for d=1:dim
        np = div(np, d)
    end
    return :($np)
end

function makeLatticeLogic(a, b, depth)
    if depth == 0
        return []
    elseif depth == 1
        results = []
        for ii = a:b
            push!(results, [ii])
        end
        return results
    else
        results = []
        for ii = a:b
            smaller_guys = makeLatticeLogic(a, b-ii+1, depth-1)
            for jj in smaller_guys
                push!(results, vcat([ii], jj))
            end
        end
        return results
    end
end

@generated function points{S<:Simplex, T<:Val}(typ::Type{EquispacedLattice{S,T}})
    d = spatialDimension(S)
    vs = getVertexCoords(S)
    n = val2val(T)
    hs = [(vs[:,i+1] - vs[:,1]) / n for i=1:d]
    icoords = makeLatticeLogic(1, n+1, d)
    results = zeros(d, size(type2type(typ)))
    pt_cur = 1
    for ii in icoords
        coord_cur = vs[:,1]
        for i=1:length(ii)
            coord_cur = coord_cur + (ii[i]-1) * hs[d+1-i]
        end
        results[:,pt_cur] = coord_cur
        pt_cur += 1
    end
    
    return :($results)
end

# tensor products
@generated function size{P1<:PointSet, P2<:PointSet}(::Type{TensorProductPointSet{P1,P2}}) 
    s = size(P1) * size(P2)
    return :($s)
end

@generated function cellType{P1<:PointSet, P2<:PointSet}(::Type{TensorProductPointSet{P1,P2}})
    ct = cellType(P1) * cellType(P2)
    return :($ct)
end

@generated function points{P1<:PointSet, P2<:PointSet}(::Type{TensorProductPointSet{P1,P2}})

end

end
