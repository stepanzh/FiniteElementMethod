#=
Reference 
    Mats G. Larson, Fredrik Bengzon.
    The Finite Element Method: Theory, Implementation, and Applications.
    Springer, 2013.
    ISBN 978-3-642-33286-9
    
    Chapter 1
=#

using DelimitedFiles


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

"Mass matrix (including boundaries) from mesh `x` and hat basis."
function mass_assemlber(x::AbstractVector{T}) where {T}
    n = length(x) - 1  # number of subintervals
    M = zeros(T, (n+1, n+1))
    @inbounds for i in 1:n  # loop over subintervals
        h = x[i+1] - x[i]
        M[i, i] += h / 3
        M[i, i+1] += h / 6
        M[i+1, i] += h / 6
        M[i+1, i+1] += h / 3
    end
    return M
end

"Trapezoidal quadrature is used."
function load_assembler(x::AbstractVector{T}, f::Function) where {T}
    n = length(x) - 1  # number of subintervals
    b = zeros(n + 1)
    @inbounds for i in 1:n
        h_half = (x[i+1] - x[i]) / 2
        b[i] += h_half * (f(x[i]))
        b[i+1] += h_half * (f(x[i+1]))
    end
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

"Dumps to `io` `xpoints`, p_h(x), p_true(x) and absolute error."
function dump(;
        io::IO=stdout,
        xpoints=collect(0:0.01:1),
        solution_fem::FEMRepresentation,
        solution_true::Function
    )
    println(io, "# x\tp_h(x)\tp_true(x)\tΔP")
    for x in xpoints
        P_fem = solution_fem(x)
        P_true = solution_true(x)
        ΔP = abs(P_true - P_fem)
        println(io, join([x, P_fem, P_true, ΔP], '\t'))
    end
    return nothing
end

function l2_projection(; mesh_size=10, target_function::Function, domain=(0, 1))
    "Uniform mesh."
    x = let n = mesh_size, Ω = domain
        h = (Ω[2] - Ω[1]) / n
        collect(Ω[1]:h:Ω[2])
    end

    M = mass_assemlber(x)
    b = load_assembler(x, target_function)
    ξ = M \ b

    "Hat basis."
    φ = basis_assembler(x)
    "Represenation."
    l²projection = FEMRepresentation(φ, ξ)

    return l²projection
end

function error_asymptotic(; target_function::Function, domain=(0, 1), io=stdout)
    mesh_error = let n = 101, Ω = domain
        h = (Ω[2] - Ω[1]) / n
        collect(Ω[1]:h:Ω[2])
    end
    
    target_discrete = map(target_function, mesh_error)
    for n in (11, 51, 101, 201, 401, 601, 801, 1001, 1201)
        solution = l2_projection(mesh_size=n, target_function=target_function, domain=domain)
        fem_discrete = map(solution, mesh_error)
        Δ = maximum(abs.(target_discrete .- fem_discrete))
        println(join([n, Δ], '\t'))
    end
end


Ω = (0, π)
target_function(x) = x^2 * sin(20x)
solution = l2_projection(
    mesh_size = 101,
    target_function = target_function,
    domain = Ω
)
dump(solution_fem=solution, solution_true=target_function, xpoints=collect(Ω[1]:1e-2:Ω[2]))    
# error_asymptotic(target_function=target_function, domain=Ω)
