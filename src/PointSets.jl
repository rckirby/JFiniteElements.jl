module PointSets

using Cells

export EquispacedLattice, Lattice, PointSet, latticeSize

abstract PointSet

abstract Lattice{C<:Cell, S<:Val} <: PointSet

immutable EquispacedLattice{C<:Cell, S<:Val} <: Lattice end

function latticeSize{C<:Cell,S<:Val}(::Type{EquispacedLattice{C,S}})
    dim = getSpatialDimension(C)
    n = val2val(S)
    return latticeSize(dim, n)
end

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
