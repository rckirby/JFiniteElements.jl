using JFiniteElements
using JFiniteElements: makepoissontri, TriangleMidpointRule


M = boxMesh(0, 1, 0, 1, 100, 100)
Q = TriangleMidpointRule()
f = makepoissontri(Q)

function buildit(M, f)
    nT = size(M.T, 2)
    cellcoords = zeros(2, 3)
    Aels = zeros(3,3,nT)
    rows = zeros(Int, 3, 3, nT)
    cols = zeros(Int, 3, 3, nT)
    TT = M.T

    for c=1:nT
        # extract/pack cell vertex coordinates
        for i=1:3
            for j=1:2
                @inbounds cellcoords[j,i] = M.x[j, M.T[i,c]]
            end
        end
        # local stiffness matrix
        f(cellcoords, sub(Aels,:,:,c))
        # and global/local coords
        for i=1:3
            for j=1:3
                @inbounds rows[i,j,c] = TT[i,c]
                @inbounds cols[i,j,c] = TT[j,c]
            end
        end
    end
    
    A = sparse(reshape(rows,3*3*nT), reshape(cols,3*3*nT), reshape(Aels, 3*3*nT))
    return A
end

function buildit2(M, f)
    nT = size(M.T, 2)
    cellcoords = zeros(2, 3)
    Ael = zeros(3,3)
#    rows = zeros(Int, 3, 3, nT)
#    cols = zeros(Int, 3, 3, nT)
    TT = M.T
    A = spzeros(size(M.x, 2), size(M.x, 2))

    for c=1:nT
        # extract/pack cell vertex coordinates
        for i=1:3
            for j=1:2
                cellcoords[j,i] = M.x[j, M.T[i,c]]
            end
        end
        # local stiffness matrix
        for i=1:3
            for j=1:3
                Ael[i,j] = 0.0
            end
        end
        f(cellcoords, Ael)
        # and global/local coords
        for i=1:3
            for j=1:3
                A[M.T[i,c], M.T[j,c]] += Ael[i,j]
            end
        end
    end
    
    return A
end

function buildit3(M, f)
    nT = size(M.T, 2)
    cellcoords = zeros(2, 3)
    Ael = zeros(3,3)
    TT = M.T


    for c=1:nT
        # extract/pack cell vertex coordinates
        for i=1:3
            for j=1:2
                @inbounds cellcoords[j,i] = M.x[j, M.T[i,c]]
            end
        end
        # local stiffness matrix
        for i=1:3
            for j=1:3
                @inbounds Ael[i,j] = 0.0
            end
        end
        f(cellcoords, Ael)
        for i=1:3
            Irow = TT[i,c]
            for j=1:3
                Jcol = TT[j,c]
                if haskey(A[Irow], Jcol)
                    @inbounds A[Irow][Jcol] += Ael[i,j]
                else
                    @inbounds A[Irow][Jcol] = Ael[i,j]
                end
            end
        end
    end
    return A
end

@time buildit(M, f)
#@time buildit2(M, f)
@time buildit3(M, f)
