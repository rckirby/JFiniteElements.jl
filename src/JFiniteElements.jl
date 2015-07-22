module JFiniteElements

export Triangulation, boxMesh

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

end # module
