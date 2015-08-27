using Cells, PointSets

C = UFCSimplex{Val{3}}
@assert getSpatialDimension(C) == 3

E = EquispacedLattice{C, Val{2}}
@assert latticeSize(E) == 10
