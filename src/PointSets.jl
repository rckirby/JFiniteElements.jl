module PointSets

using Cells

export EquispacedLattice, Lattice, PointSet, latticeSize

abstract PointSet

abstract Lattice{C<:Cell, S<:Val} <: PointSet

immutable EquispacedLattice{C<:Cell, S<:Val} <: Lattice end

rootleaf{S,T}(::Type{EquispacedLattice{S,T}}) = Lattice{S,T}


latticeSize{T}(::Type{T}) = latticeSize(rootleaf(T))

# This can be "generated" since it just needs to be run once per type.
# After that, the first-time value will be returned.
@generated function latticeSize{S<:Simplex,T<:Val}(::Type{Lattice{S,T}})
    println("hi there")
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

# Doesn't work yet.  Need an iterator over lattices in place first.
function latticePoints{S<:Simplex, T<:Val}(::Type{EquispacedLattice{S,T}})
    d = getSpatialDimension(S)
    vs = getVertexCoords(S)
    n = val2val(T)
    hs = [(vs[i+1] - vs[1]) / n for i in range(1,d+1)]
    
    return
end

end
