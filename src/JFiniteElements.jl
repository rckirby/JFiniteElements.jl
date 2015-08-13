module JFiniteElements

using Bases, Cells
export Bases, Cells

export Triangulation, boxMesh, poissontri

type Triangulation
    x    
    T
end

function boxMesh(a, b, c, d, Nx, Ny)
    xvals = linspace(a, b, Nx+1)
    yvals = linspace(c, d, Ny+1)

    x::Array{Float64, 2} = zeros(2, (Nx+1)*(Ny+1))
    T::Array{Int, 2} = zeros(3, 2*Nx*Ny)

    for j=1:Nx+1
        for i=1:Ny+1
            idx = (j-1)*(Ny+1) + i
            x[1, idx] = xvals[j]
            x[2, idx] = yvals[i]
        end
    end

    for j=1:Nx
        for i=1:Ny
            # Global square index
            S = i + (j-1)*Ny
            # Two global triangles coming from square S
            T1 = 2S - 1
            T2 = 2S
            # Four vertices of the square
            # numbered by SW, NW, SE, NE
            v1 = i+(j-1)*(Ny+1)
            v2 = v1 + 1
            v3 = v1 + Ny+1
            v4 = v3 + 1
            # Three vertices of each triangle
            T[1, T1] = v1
            T[2, T1] = v2
            T[3, T1] = v3
            T[1, T2] = v2
            T[2, T2] = v3
            T[3, T2] = v4
        end
    end

    return Triangulation(x, T)
end

#######
# Experimental type system.  We'll see if it works...
#
abstract Cell

abstract Simplex <: Cell

# We could get fancy and Vararg this, but pairs work fine?
type TensorProduct{A, B} <: Cell
    Cell1::A
    Cell2::B
end

abstract CellMapping

abstract VertexMapping <: CellMapping

abstract AffineMapping <: VertexMapping

abstract PointSet

# Note, I *don't* want to just have a single type
# since then I can't multiply dispatch evaluation algorithms
# at points.  I could parameterize the type?
abstract QuadratureRule

immutable TriangleMidpointRule <: QuadratureRule end
immutable TriangleVertexRule <: QuadratureRule end
immutable TriangleEdgeRule <: QuadratureRule end

getQuadratureWtsPts(T::TriangleMidpointRule) = 
  ( [0.5], reshape([1./3, 1./3], 2, 1))

numQuadraturePoints(T::TriangleMidpointRule) = 1

getQuadratureWtsPts(T::TriangleVertexRule) =
    ( [1./6, 1./6, 1./6], reshape([0 0 1 0 0 1.], 2, 3))

numQuadraturePoints(T::TriangleVertexRule) = 3

getQuadratureWtsPts(T::TriangleVertexRule) =
    ([1./6, 1./6, 1./6], reshape([0 0 1 0 0 1.], 2, 3))

numQuadraturePoints(T::TriangleEdgeRule) = 3

getQuadratureWtsPts(T::TriangleEdgeRule) =
    ([1./6, 1./6, 1./6], reshape([0.5 0 0.5 0.5 0 0.5], 2, 3))


#
# End experimental type system.
#######



###
# Let's test out and refine the type system by considering some special
# use cases we'd like
#


# This is Poisson with P1 on a particular reference element
# with the midpoint rule for quadrature.

function poissontri(x, A)
    Nqp = 1
    Nbf = 3

    # First, compute the Jacobian.  In this case, it's constant
    J = reshape([[x[1,2] - x[1,1], x[2,2] - x[2,1]];
                 [x[1,3] - x[1,1], x[2,3] - x[2,1]]], 2, 2)
    detJ = det(J)
    Jinv = inv(J)

    # Next, tabulate the basis function gradients at qp.  We'll have
    # gradients are two-dimensional here, there is one quadrature point,
    # and there are three basis functions.  That explains the size of
    # this array 
    DPsi = zeros(2, Nqp, Nbf)
    DPsi[:,1,1] = [-1; -1]
    DPsi[:,1,2] = [1; 0]
    DPsi[:,1,3] = [0; 1]

    DPsiPhysQP = zeros(2, Nbf)

    w = [0.5]
    
    for q=1:Nqp
        # transform gradients at q:th quadrature point
        for i=1:Nbf
            DPsiPhysQP[:,i] = Jinv' * DPsi[:,q,i]
        end
        for i=1:Nbf
            for j=1:Nbf
                A[i, j] += w[q] * dot(DPsiPhysQP[:,i], DPsiPhysQP[:, j])
            end
        end
    end
    # Since the Jacobian determinant is constant, we apply it here.
    A[:, :] *= abs(detJ)
end

function poissonquad(x, A)
    Nqp = 1
    Nbf = 4

    # Note that the Jacobians are spatially varying now, so we need
    # the quadrature point locations on the reference cell.  Still
    # using the midpoint rule
    quadPts = reshape([0.5 0.5], 2, 1)
    quadWts = [1.0]

    # unpack to make programming easier
    x1 = x[1,1]
    y1 = x[2,1]
    x2 = x[1,2]
    y2 = x[2,2]
    x3 = x[1,3]
    y3 = x[2,3]
    x4 = x[1,4]
    y4 = x[2,4]
    x21 = x2-x1
    y21 = y2-y1
    x31 = x3-x1
    y31 = y3-y1
    x4123 = x4+x1-(x2+x3)
    y4123 = y4+y1-(y2+y3)
    
    J = zeros(2, 2, Nqp)
    Jinv = zeros(2, 2, Nqp)
    detJ = zeros(Nqp)
    
    for q=1:Nqp
        xhat = quadPts[1,q]
        yhat = quadPts[2,q]
        J[1, 1, q] = x21 + x4123*yhat
        J[2, 1, q] = y21 + y4123*yhat
        J[1, 2, q] = x31 + x4123*xhat
        J[2, 2, q] = y31 + y4123*xhat
        Jinv[:, :, q] = inv(J[:, :, q])
        detJ[q] = det(J[:, :, q])
    end
    # Done with Jacobian.

    # Next, tabulate the basis function gradients at qp.  We'll have
    # gradients are two-dimensional here, there is one quadrature point,
    # and there are four basis functions.  That explains the size of
    # this array
    # BF are:
    # (1-xhat)(1-yhat)
    # xhat(1-yhat)
    # yhat(1-xhat)
    # xhat yhat
    # so that gradidents are
    # [yhat-1; xhat-1]
    # [1-yhat; -xhat]
    # [-yhat; 1-xhat]
    # [yhat; xhat]
    DPsi = zeros(2, Nqp, Nbf)
    for q=1:Nqp
        xhat = quadPts[1, q]
        yhat = quadPts[2, q]
        DPsi[1, q, 1] = yhat-1
        DPsi[2, q, 1] = xhat-1
        DPsi[1, q, 2] = 1-yhat
        DPsi[2, q, 2] = -xhat
        DPsi[1, q, 3] = -yhat
        DPsi[2, q, 3] = 1-xhat
        DPsi[1, q, 4] = yhat
        DPsi[2, q, 4] = xhat
    end
    

    DPsiPhysQP = zeros(2, Nbf)

    
    for q=1:Nqp
        # transform gradients at q:th quadrature point
        for i=1:Nbf
            DPsiPhysQP[:,i] = Jinv[:,:,q]' * DPsi[:,q,i]
        end

        for i=1:Nbf
            for j=1:Nbf
                A[i, j] += abs(detJ[q]) * quadWts[q] * dot(DPsiPhysQP[:,i], DPsiPhysQP[:, j])
            end
        end
    end
end

# This fixes some of the allocation problems by not using the operators
# and writing explicit loops for the transformation/dot products

function poissontri2(x, A)
    Nqp = 1
    Nbf = 3

    # First, compute the Jacobian.  In this case, it's constant
    J = reshape([[x[1,2] - x[1,1], x[2,2] - x[2,1]];
                 [x[1,3] - x[1,1], x[2,3] - x[2,1]]], 2, 2)
    detJ = det(J)
    Jinv = inv(J)

    # Next, tabulate the basis function gradients at qp.  We'll have
    # gradients are two-dimensional here, there is one quadrature point,
    # and there are three basis functions.  That explains the size of
    # this array 
    DPsi = zeros(2, Nqp, Nbf)
    DPsi[:,1,1] = [-1; -1]
    DPsi[:,1,2] = [1; 0]
    DPsi[:,1,3] = [0; 1]

    DPsiPhysQP = zeros(2, Nbf)

    w = [0.5]
    
    for q=1:Nqp
        # transform gradients at q:th quadrature point
        for i=1:Nbf
            for j=1:2
                DPsiPhysQP[j, i] = 0.0
                for k=1:2
                    DPsiPhysQP[j,i] += Jinv[k, j] * DPsi[k,q,i]
                end
            end
        end
        for i=1:Nbf
            for j=1:Nbf
                for k=1:2
                    A[i, j] += w[q] * DPsiPhysQP[k,i] * DPsiPhysQP[k, j]
                end
            end
        end
    end
    absdetJ = abs(detJ)
    for j=1:Nbf
        for i=1:Nbf
            A[i,j] *= absdetJ
        end
    end
end

#
########

# This is a stab at pre-setting the quadrature rule, basis function values, and preallocating space in Poisson on a triangle.

function makepoissontri(Q)
    Nqp = numQuadraturePoints(Q)
    w, quadpts = getQuadratureWtsPts(Q)
    Nbf = 3

    # Allocate space for Jacobian and its inverse and determinant

    J = zeros(2, 2)
    adetJ::Float64 = 0.0
    detJ::Float64 = 0.0
    Jinv = zeros(2, 2)
    i::Int = 0
    j::Int = 0
    k::Int = 0
    q::Int = 0

    # quadrature information is baked in.  One-point midpoint rule

    # Next, tabulate the basis function gradients at qp.  We'll have
    # gradients are two-dimensional here, there is one quadrature point,
    # and there are three basis functions.  That explains the size of
    # this array 
    DPsi = zeros(2, Nqp, Nbf)
    for q=1:Nqp
        @inbounds DPsi[:,q,1] = [-1; -1]
        @inbounds DPsi[:,q,2] = [1; 0]
        @inbounds DPsi[:,q,3] = [0; 1]
    end

    DPsiPhysQP = zeros(2, Nbf)

    function f(x, A)
        # Write the Jacobian into pre-allocated space
        @inbounds J[1,1] = x[1,2] - x[1,1]
        @inbounds J[1,2] = x[2,2] - x[2,1]
        @inbounds J[2,1] = x[1,3] - x[1,1]
        @inbounds J[2,2] = x[2,3] - x[2,1]
        @inbounds detJ = J[1,1] * J[2,2] - J[1,2] * J[2,1]
        @inbounds Jinv[1,1] = J[2,2] / detJ
        @inbounds Jinv[1,2] = -J[1,2] / detJ
        @inbounds Jinv[2,1] = -J[2,1] / detJ
        @inbounds Jinv[2,2] = J[1,1] / detJ

        for q=1:Nqp
            # transform gradients at q:th quadrature point
            for i=1:Nbf
                for j=1:2
                    DPsiPhysQP[j, i] = 0.0
                    for k=1:2
                        @inbounds DPsiPhysQP[j,i] += Jinv[k, j] * DPsi[k,q,i]
                    end
                end
            end
            for i=1:Nbf
                for j=1:Nbf
                    for k=1:2
                        @inbounds A[i, j] += w[q] * DPsiPhysQP[k,i] * DPsiPhysQP[k, j]
                    end
                end
            end
        end
        # Since the Jacobian determinant is constant, we apply it here.
        adetJ = detJ < 0 ? -detJ : detJ
        for i=1:Nbf
            for j=1:Nbf
               @inbounds A[i, j] *= adetJ
            end
        end
    end

    
    return f
end


end # module
