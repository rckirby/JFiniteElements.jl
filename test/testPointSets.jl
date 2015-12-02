using JFiniteElements, Cells, PointSets

T = Cells.UFCSimplex{Val{2}}
EL = PointSets.EquispacedLattice{T,Val{3}}
@assert PointSets.latticeSize(EL) == 10

