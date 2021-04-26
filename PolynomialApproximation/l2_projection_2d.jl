#=
Reference 
    Mats G. Larson, Fredrik Bengzon.
    The Finite Element Method: Theory, Implementation, and Applications.
    Springer, 2013.
    ISBN 978-3-642-33286-9
    
    Chapter 3.6
=#

using SparseArrays, LinearAlgebra

"Stores point, connectivity matrices and satellite data."
struct Triangulation{F,I}
    P::Matrix{F}
    T::Matrix{I}
    nxy::Vector{I}
    hxy::Vector{F}
end

point_matrix(x::Triangulation) = x.P
connectivity_matrix(x::Triangulation) = x.T
number_points(x::Triangulation) = size(point_matrix(x), 2)
number_elements(x::Triangulation) = size(connectivity_matrix(x), 2)

"Returns coordinates of vertexes of element `K`."
function element_vertexes(mesh::Triangulation, K::Int)
    N = connectivity_matrix(mesh)[:, K]
    vertexes = [point_matrix(mesh)[:, K] for K in N]
    return vertexes
end

"Returns area of triangle `K`."
function area(mesh::Triangulation, K::Int)
    "Vertex indices."
    vertexes = element_vertexes(mesh, K)
    AB = vertexes[1] - vertexes[2]
    BC = vertexes[2] - vertexes[3]
    S = 0.5 * abs(AB[1] * BC[2] - AB[2]*BC[1])  # cross product in 2D used
    return S
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

"Returns mass matrix for a basis of 2d hat-functions."
function mass_assembler_2d(mesh::Triangulation{F,I}) where {F, I}
    P, T = point_matrix(mesh), connectivity_matrix(mesh)
    np, nt = number_points(mesh), number_elements(mesh)
    M = spzeros(F, np, np)
    for K in 1:nt
        loc2glb = T[1:3, K]
        K_area = area(mesh, K)
        MK = [2 1 1; 1 2 1; 1 1 2] ./ 12 .* K_area
        # M[r, r] += MK[1, 1], M[r, s] += MK[1, 2], ...
        M[loc2glb, loc2glb] += MK
    end
    return M
end

"""
Returns load vector for 2d hat-functions basis and function `f`.
Uses node quadrature rule.
"""
function load_assembler_2d(mesh::Triangulation{F,I}, f::Function) where {F,I}
    b = zeros(F, number_points(mesh))
    for K in 1:number_elements(mesh)
        loc2glb = connectivity_matrix(mesh)[1:3, K]
        K_area = area(mesh, K)
        N = element_vertexes(mesh, K)
        bK = 1/3 .* map(f, N) .* K_area  # node-quadrature
        b[loc2glb] += bK
    end
    return b
end


target_function(x) = x[1]^2 + x[2]^2

mesh = TriangulationOfRect([-1, -1], [1, 1], [10, 10])
M = mass_assembler_2d(mesh)
b = load_assembler_2d(mesh, target_function)
ξ = M \ b

for i in eachindex(ξ)
    Nᵢ = point_matrix(mesh)[:, i]
    πf = ξ[i]
    f = target_function(Nᵢ)
    println(join([Nᵢ..., f, πf, f - πf], '\t'))
end
