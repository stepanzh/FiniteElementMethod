#=
Reference 
    Mats G. Larson, Fredrik Bengzon.
    The Finite Element Method: Theory, Implementation, and Applications.
    Springer, 2013.
    ISBN 978-3-642-33286-9
    
    Chapters 2.1 - 2.4

Problem

  Domain Ω = [0, L]
  ODE
    -(a(x) u')' = f(x), x ∈ (0, L), a(x) > 0

  Boundary conditions
    a(0)u'(0) = κ₀ (u(0) - g₀)
    -a(L)u'(L) = κ₁ (u(L) - g₁)

   Parameters
    κ₀ ≥ 0, κ₁, gₒ, g₁
=#

"Hat basis function φⱼ."
struct Hat{T<:Real} <: Function
    x₋::T
    x₀::T
    x₊::T
end
function (f::Hat{T})(x::Real) where {T}
    if x == f.x₀  # is ok to use `==` here (?)
        return one(T)
    elseif f.x₋ < x < f.x₀
        Δx = f.x₀ - f.x₋
        return (x - f.x₋) / Δx
    elseif f.x₀ < x < f.x₊
        Δx = f.x₊ - f.x₀
        return (f.x₊ - x) / Δx
    else
        return zero(T)
    end
end

"υ(x) = ∑υᵢφᵢ(x)"
struct FEMRepresentation{T}
    φ::Vector{Hat{T}}
    υ::Vector{T}
end
function (f::FEMRepresentation{T})(x::Real) where {T}
    val = zero(T)
    for (φᵢ, υᵢ) in zip(f.φ, f.υ)
        val += υᵢ * φᵢ(x)
    end
    return val
end

function stiffness_assembler(x, conductivity, κ∂Ω)
    n = length(x) - 1
    A = zeros((n + 1, n + 1))
    @inbounds for i in 1:n
        h = x[i+1] - x[i]
        xmid = (x[i+1] + x[i]) / 2
        amid = conductivity(xmid)

        A[i, i] += amid / h
        A[i, i+1] -= amid / h
        A[i+1, i] -= amid / h
        A[i+1, i+1] += amid / h
    end
    A[1, 1] += κ∂Ω[1]
    A[n+1, n+1] += κ∂Ω[2]
    return A
end

function source_assembler(x, heatsource, κ∂Ω, g∂Ω)
    n = length(x) - 1
    b = zeros(n + 1)
    @inbounds for i in 1:n
        h_half = (x[i+1] - x[i]) / 2
        b[i] += h_half * (heatsource(x[i]))
        b[i+1] += h_half * (heatsource(x[i+1]))
    end
    b[1] += κ∂Ω[1] * g∂Ω[1]
    b[n+1] += κ∂Ω[2] * g∂Ω[2]
    return b
end

function basis_assembler(x::AbstractVector{T}) where {T}
    φ = Vector{Hat{T}}(undef, length(x))
    φ[1] = Hat(NaN, x[1], x[2])  # half-hat
    @inbounds for i in 2:length(x)-1
        φ[i] = Hat(x[i-1], x[i], x[i+1])
    end
    φ[end] = Hat(x[end-1], x[end], NaN)  # half-hat
    return φ
end

"Uniform mesh of Ω including Ω[1] and Ω[2] at x[1] and x[n] respectively."
function mesh_assembler(n::Int, Ω)
    h = (Ω[2] - Ω[1]) / (n - 1)
    return collect(Ω[1]:h:Ω[2])
end

"""
Solver of stationary parabolic problem.

- `Ω`: domain, must behave as tuple of two elements (0, L)
- `κ∂Ω`: boundary κ, must behave as tuple (κ(0), κ(L))
"""
function solver(; mesh_points::Int=25, Ω, κ∂Ω, g∂Ω, conductivity::Function, heatsource::Function)
    "Mesh (partition) of Ω."
    x = mesh_assembler(mesh_points, Ω)

    stiffness = stiffness_assembler(x, conductivity, κ∂Ω)
    source_vector = source_assembler(x, heatsource, κ∂Ω, g∂Ω)
    uᵢ = stiffness \ source_vector
    φ = basis_assembler(x)
    uh = FEMRepresentation(φ, uᵢ)
    return uh
end

function dump(; io::IO=stdout, mesh_points=25, Ω, solution::FEMRepresentation)
    mesh = mesh_assembler(mesh_points, Ω)
    println(io, "# x\tsolution(x)")
    for x in mesh
        println(io, join([x, solution(x)], '\t'))
    end
end

# Problem statement
Ω = (2, 8)
κ = (1e6, 0)
@assert any(x -> x > 0, κ) "negative κ"
g = (-1, 0)

"Function a(x) > 0, x ∈ Ω."
conductivity(x::Real) = 0.1 * (5 - 0.6x)
"Function f(x), x ∈ Ω."
heatsource(x::Real) = 0.03 * (x - 6)^4

solution_fem = solver(
    mesh_points=50,
    Ω=Ω,
    κ∂Ω=κ,
    g∂Ω=g,
    conductivity=conductivity,
    heatsource=heatsource
)
dump(
    io=stdout,
    mesh_points=100,
    Ω=Ω,
    solution=solution_fem
)
