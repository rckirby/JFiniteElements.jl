module PointSets

using Cells

export EquispacedLattice, Lattice, PointSet, latticeSize

abstract PointSet

abstract Lattice{C<:Cell, S<:Val} <: PointSet

immutable EquispacedLattice{C<:Cell, S<:Val} <: Lattice end

rootleaf{S,T}(::Type{EquispacedLattice{S,T}}) = Lattice{S,T}

# as we have other kinds of cells, we need to specify their rootleafs
# and then dispatch on sizes.
# unless we can get away with just adding via tensor products.

function latticeSize{S<:Simplex,T<:Val}(::Type{Lattice{S,T}})
    return latticeSize(getSpatialDimension(S),
                       val2val(T))
end

latticeSize{T}(::Type{T}) = latticeSize(rootleaf(T))


function latticeSize(dim, n)
    np = 1
    for d=1:dim
        np *= (n+d)
    end
    for d=1:dim
        np = div(np, d)
    end
    return np
end

end
