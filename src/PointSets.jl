module PointSets

using Cells

export EquispacedLattice, Lattice, PointSet, latticeSize

abstract PointSet

abstract Lattice{C<:Cell, S<:Val} <: PointSet

immutable EquispacedLattice{C<:Cell, S<:Val} <: Lattice end

rootleaf{S,T}(::Type{EquispacedLattice{S,T}}) = Lattice{S,T}


latticeSize{T}(::Type{T}) = latticeSize(rootleaf(T))


type2type{T}(::Type{Type{T}}) = T

# This can be "generated" since it just needs to be run once per type.
# After that, the first-time value will be returned.
@generated function latticeSize{S<:Simplex,T<:Val}(::Type{Lattice{S,T}})
    println("lattice size")
    dim = getSpatialDimension(S)
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

@generated function latticePoints{S<:Simplex, T<:Val}(typ::Type{EquispacedLattice{S,T}})
    d = getSpatialDimension(S)
    vs = getVertexCoords(S)
    n = val2val(T)
    hs = [(vs[:,i+1] - vs[:,1]) / n for i=1:d]
    icoords = makeLatticeLogic(1, n+1, d)
    results = zeros(d, latticeSize(type2type(typ)))
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

end
