M = boxMesh(0, 1, 0, 1, 3, 5)
@test length(M.x) == 2 * 24
@test length(M.T) == 3 * 2 * 15 
