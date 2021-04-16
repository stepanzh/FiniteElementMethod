#=
Reference
  Chen, Z., Huan, G., & Ma, Y. (2006). Computational methods for multiphase flows in porous media. Society for Industrial and Applied Mathematics.
  Chapter 4.2
  ISBN 0-89871-606-3

Problem
  Stationary linear elliptic equation (Poisson's)
  Domain Ω = [0, 1] (bounded)

    d²P
  - --- = f(x),   0 < x < 1
    dx²
  
  where f(x) = Cf = const

  Boundary conditions (BC)
    P(0) = P₀ = 0
    P(1) = P₁ = 0

Analytic solution of the problem (BC not only zeros)
  P(x) = - 0.5 Cf x² + C₁ x + C₂   where C₁ = P₁ + 0.5 Cf - C₂  and  C₂ = P₀

Mesh
  Partition - defined by number M > 0
    0 = x₀ < x₁ < ... < x_M < x_{M+1} = 1
  Subinterval Iᵢ
    Iᵢ = (xᵢ₋₁, xᵢ)   length(Iᵢ) = xᵢ - xᵢ₋₁   i = 1, 2, ..., M + 1

Finite Element Space V_h
  V_h = {
     function υ:
        υ is continious on Ω,
        υ is linear on each Iᵢ,
        υ(0) = υ(1) = 0
  }

Basis functions φᵢ ∈ V_h, i = 1, 2, ..., M
  φᵢ - hat (chapeau) function with top at xᵢ
  φᵢ(xᵢ) = 1
  
  1      ╱╲        x₋ = xᵢ₋₁
        ╱  ╲       x₀ = xᵢ
       ╱    ╲      x₊ = xᵢ₊₁
  0 ▔▔*▔▔*▔▔▔*▔▔
      x₋ x₀  x₊

Representation of solution P_h(x)

  P_h(x) = Σ pᵢ φᵢ(x), where pᵢ = p_h(xᵢ), the sum goes over i = 1, 2, ..., M

Values pᵢ are found from the following system

 Ap = F, where
   A - MxM 'stiffness' matrix, aᵢⱼ = (φᵢ, φⱼ)
   p - Mx1 vector of pᵢ
   F - Mx1 'source' vector of (f, φᵢ)
=#


"Subinterval Iᵢ."
struct Subinterval{T<:Real}
    x₁::T
    x₂::T
end
left(x::Subinterval) = x.x₁
right(x::Subinterval) = x.x₂
len(x::Subinterval) = right(x) - left(x)

"Hat basis function φⱼ."
struct Hat{T<:Real}
    x₋::T
    x₀::T
    x₊::T
end
fvalue(f::Hat{T}, x::Real) where {T} = begin
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
    return nothing  # for errors
end

"υ(x) = ∑υᵢφᵢ(x)"
struct FEMRepresentation{T}
    φ::Vector{Hat{T}}
    υ::Vector{T}
end
fvalue(f::FEMRepresentation{T}, x::Real) where {T} = begin
    val = zero(T)
    for (φᵢ, υᵢ) in zip(f.φ, f.υ)
        val += fvalue(φᵢ, x) * υᵢ
    end
    return val
end

"Stiffness matrix for hat basis functions φᵢ."
function stiffness_matrix(S::AbstractVector{Subinterval{T}}) where {T}
    M = length(S) - 1
    A = zeros((M, M))
    
    A[1, 1] = 1 / len(S[1]) + 1 / len(S[2])
    A[1, 2] = - 1 / len(S[2])
    
    A[M, M-1] = - 1 / len(S[M])
    A[M, M] = 1 / len(S[M]) + 1 / len(S[M+1])

    @inbounds for i in 2:size(A)[1]-1
          hᵢ = len(S[i])
        hᵢ₊₁ = len(S[i+1])
        A[i, i-1] = - 1 / hᵢ
        A[i, i] = (1 / hᵢ) + (1 / hᵢ₊₁)
        A[i, i+1] = - 1 / hᵢ₊₁
    end
    return A
end

"Source vector Fⱼ = (f(x), φⱼ), j = 1, ..., M."
function source_vector_hat_constant(S::AbstractVector{Subinterval{T}}, Cf::Real=1.0) where {T}
    F = fill(NaN, length(S) - 1)
    for i in eachindex(F)
        F[i] = 0.5 * Cf * (len(S[i]) + len(S[i+1]))
    end
    return F
end

"Generates Vector of hat functions (basis)."
function generate_hat_basis(S::AbstractVector{Subinterval{T}}) where {T}
    φ = Vector{Hat{T}}(undef, length(S)-1)
    @inbounds for i in eachindex(φ)
        I₋ = S[i]
        I₊ = S[i+1]
        x₋, x₀, x₊ = left(I₋), right(I₋), right(I₊)  # right(I₋) == left(I₊)
        φ[i] = Hat(x₋, x₀, x₊)
    end
    return φ
end

"Analytic (true) solution of the problem. Returns closure."
function analytic_solution(P₀::Real, P₁::Real, Cf::Real)
    C₂ = P₀
    C₁ = P₁ + Cf / 2 - C₂
    return x -> - 0.5 * Cf * x^2 + C₁ * x + C₂
end

"Dumps to `io` `xpoints`, p_h(x), p_true(x) and absolute error."
function dump(;
        io::IO=stdout,
        xpoints=collect(0:0.01:1),
        solution_fem::FEMRepresentation,
        solution_true::Function
    )
    println(io, "# x\tp_h(x)\tp_true(x)\tΔP")
    for x in xpoints
        P_fem = fvalue(solution_fem, x)
        P_true = solution_true(x)
        ΔP = abs(P_true - P_fem)
        println(io, join([x, P_fem, P_true, ΔP], '\t'))
    end
    return nothing
end

"Generates mesh based on its `type` and number of inner points `M`."
function generate_mesh(; M::Integer, type="linear")
    if type == "linear"
        h = 1 / (M + 1)
        return collect(0.0:h:1.0)
    elseif type == "random"
        mesh = Vector{Float64}(undef, M+2)
        mesh[1] = 0
        mesh[2:M+1] .= sort(rand(M))
        mesh[M+2] = 1
        return mesh
    else
        @error "Unrecognized `type` of mesh."
    end
end

# Problem statement
"Boundary."
P₀ = P₁ = 0
"f(x) = Cf."
Cf = 1

# Partition of [0, 1]. xᵢ points i=0, 1, ..., M+1
"xᵢ points i = 0, 1, ..., M+1"
mesh = generate_mesh(M=199, type="random")

"Vector of subintervals Iᵢ."
intervals = let mesh=mesh
    map(i -> Subinterval(mesh[i], mesh[i+1]), 1:length(mesh)-1)
end

"Stiffness matrix."
A = stiffness_matrix(intervals)
"Source vector."
F = source_vector_hat_constant(intervals, Cf)
"Pᵢ vector."
Pᵢ = A \ F
"Vector of basis functions."
φ_basis = generate_hat_basis(intervals)

"FEM solution."
solution_repr = FEMRepresentation(φ_basis, Pᵢ)
"Analytic solution."
solution_true = analytic_solution(P₀, P₁, Cf)

dump(
    solution_fem=solution_repr,
    solution_true=solution_true
)
