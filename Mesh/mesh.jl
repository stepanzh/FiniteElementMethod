#=
Reference 
    Mats G. Larson, Fredrik Bengzon.
    The Finite Element Method: Theory, Implementation, and Applications.
    Springer, 2013.
    ISBN 978-3-642-33286-9
    
    Chapter 3.1
=#

"Stores point, connectivity matrices and satellite data."
struct Triangulation{F,I}
    P::Matrix{F}
    T::Matrix{I}
    nxy::Vector{I}
    hxy::Vector{F}
end

"""
Returns `Triangulation` of rectangular area defined by
bottom-left corner `p₁`, top-right corner `p₂`
and partition of x and y side `nxy`.

All arguments must behave as `Vector` in 2D.
"""
function TriangulationOfRect(p₁::V, p₂::V, nxy::V) where {V <: AbstractVector}
    hxy = (p₂ .- p₁) ./ nxy
    nx, ny = nxy
    hx, hy = hxy

    #=
    Cycle goes over all (i, j) points.
    =#
    "Point matrix."
    P = Matrix{Float64}(undef, (2, (nx + 1)*(ny + 1)))
    @inbounds for i in 1:(nx + 1)
        for j in 1:(ny + 1)
            node = i + (j - 1)*(nx + 1)        # number of node from (i, j)
            P[1, node] = p₁[1] + (i - 1) * hx  # coordinates of node
            P[2, node] = p₁[2] + (j - 1) * hy  #
        end
    end

    #=
    Connectivity matrix

    (i, j) cycle goes       o─o─o─o─o
    from bottom to top,     │╲│╲│╲│╲│
    from left to right,     *─*─*─*─o
    on *-point              │╲│╲│╲│╲│
                            *─*─*─*─o
                            │╲│╲│╲│╲│
                            *─*─*─*─o

    Cycle body looks at     N₃──N₄    where N₁ is number of (i,j)-node
                element     │ ╲ │
                            N₁──N₂
    =#
    "Connectivity matrix."
    T = Matrix{Int}(undef, (3, 2*nx*ny))
    @inbounds for i in 1:nx
        for j in 1:ny
            K₁ = 2 * i + 2 * nx * (j - 1) - 1
            K₂ = K₁ + 1
            N₁ = i + (j - 1)*(nx + 1)
            N₂ = N₁ + 1
            N₃ = N₁ + nx + 1
            N₄ = N₃ + 1
            T[:, K₁] = [N₁, N₂, N₃]
            T[:, K₂] = [N₂, N₄, N₃]
        end
    end
    return Triangulation{Float64,Int}(P, T, [nx, ny], [hx, hy])
end

function dump(K::Triangulation)
    for Kⱼ in 1:size(K.T)[2]
        println("# K_$Kⱼ")
        for vertex in K.T[:, Kⱼ]
            println(join(K.P[:, vertex], '\t'))
        end
        print("\n\n\n")
    end
end

K = TriangulationOfRect([0, 0], [1, 1], [4, 4])
dump(K)
