module Bases

abstract Basis

immutable P1Triangle <: Basis end

numMembers(B::P1Triangle) = 3

function getPsi(B::P1Triangle, x, Psi)
    # x should be 2 x num points
    # Psi should be num points x 3
    @assert size(x,1) == 2
    @assert size(Psi,1) == size(x,2)
    @assert size(Psi,2) == 3

    np = size(x,2)

    for i=1:np
        Psi[i,1] = 1 - x[1,i] - x[2,i]
        Psi[i,2] = x[1,i]
        Psi[i,3] = x[2,i]
    end

    return
end

function getDPsi(B::P1Triangle, x, DPsi)
    # x should be 2 x num points
    # DPsi should be 2 x num points x 3
    @assert size(x,1) == 2
    @assert size(DPsi,1) == 2
    @assert size(DPsi,2) == size(x,2)
    @assert size(DPsi,3) == 3

    np = size(x,2)

    for i=1:np
        DPsi[1,i,1] = -1.0
        DPsi[2,i,1] = -1.0
        DPsi[1,i,2] = 1.0
        DPsi[2,i,2] = 0.0
        DPsi[1,i,3] = 0.0
        DPsi[2,i,3] = 1.0
    end

    return
end

end
