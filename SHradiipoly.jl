# solve Swift-Hohenberg via Newton

function F!(F, u, λ)
    ∂² = project(Derivative(2), space(u), space(u))
    project!(F, λ * u - (I + ∂²)^2 * u - u^3)
    return F
end

function DF!(DF, u, λ)
    ∂² = project(Derivative(2), space(u), space(u))
    add!(DF, λ * I - (I + ∂²)^2, Multiplication(-3*u^2))
    return DF
end

using RadiiPolynomial

# m = 10
# N = m-1
# λ = 1.0
# L = 2π
# ω = 1.0
# u = Sequence(CosFourier(N, ω), rand(N+1))


m = 51
N = m-1
λ = 5.0
L = 6π
ω = 1.0/3.0 # 2π/L
u = Sequence(CosFourier(N, ω), 100*randn(N+1)./(1:N+1).^3)

newton!((F, DF, u) -> (F!(F, u, λ), DF!(DF, u, λ)), u, maxiter=100)


# using GLMakie
# lines(LinRange(0, π, 101), x -> u(x); linewidth = 6)

using Plots
plot(x -> u(x),0 ,L, legend=false, title = "Plot cosine Fourier series",
        line=2,
        xlabel = "\$x\$",
        ylabel = "\$u(x)\$",)


#######################################
# Start validation
#######################################

function iF!(F, u, λ)
    ∂² = project(Derivative(2), space(u), space(u))
    project!(F, λ * u - interval.((I + ∂²)^2) * u - u^3)
    return F
end

function iDF!(DF, u, λ)
    ∂² = project(Derivative(2), space(u), space(u))
    iI = interval(Matrix(1.0I, length(u), length(u)))
    μ = LinearOperator(space(u), space(u), λ * iI) - interval.((I + ∂²)^2)
    add!(DF, μ, Multiplication(-interval(3)*u^2))
    return DF
end

# A† , An
iu = interval.(u)
iλ = interval(λ)
iν = interval(1.0)

# Function space
X = Ell1(GeometricWeight(iν))
# 
_DF_ = DF!(zeros(eltype(u), CosFourier(N, ω), CosFourier(N, ω)), u, λ)
iA = interval.(inv(_DF_))
tail_A_norm = inv(abs(iλ - interval(1 - (N+1)^2)^2))

# Y0 bound
iF = iF!(zeros(eltype(iu),CosFourier(3N, ω)), iu, iλ)
tail_F = copy(iF)
tail_F[0:N] .= interval(0)
Y = norm(iA * iF, X) + norm(tail_F, X) * tail_A_norm

# Z bound
iDF = iDF!(zeros(eltype(iu),CosFourier(3N, ω), CosFourier(N, ω)), iu, iλ)

iu²3 = interval(3)iu^2
tail_u²3 = copy(iu²3)
tail_u²3[0:N] .= interval(0)

# Z0 + Z1
iI = LinearOperator(CosFourier(3N, ω), CosFourier(N, ω), interval(Matrix(1.0I,N+1,3N+1)))
Z = max(opnorm(iA * iDF - iI, X) + norm(tail_u²3, X) * tail_A_norm,
        norm(iu²3, X) * inv(abs(iλ - interval(1 - (3N+1)^2)^2)))

#

R = interval(2sup(Y))
W = interval(6) * max(opnorm(iA, X), tail_A_norm)*(norm(iu, X) + R)
Z2 = interval(3) * max(opnorm(iA, X), tail_A_norm)*(interval(2) * norm(iu, X) + R)
#
r = interval_of_existence(Y, Z, W, mid(R))



# # Multidim

# function F!(F, u, λ)
#     ∂1² = project(Derivative(2, 0), space(u), space(u))
#     ∂2² = project(Derivative(0, 2), space(u), space(u))
#     Δ = ∂1² + ∂2²
#     project!(F, λ * u - (I + Δ)^2 * u - u^3)
#     return F
# end

# function DF!(DF, u, λ)
#     ∂1² = project(Derivative(2, 0), space(u), space(u))
#     ∂2² = project(Derivative(0, 2), space(u), space(u))
#     Δ = ∂1² + ∂2²
#     add!(DF, λ * I - (I + Δ)^2, Multiplication(-3u^2))
#     return DF
# end

# using RadiiPolynomial

# K = 20
# λ = 1.0
# u = Sequence(CosFourier(K, 0.5) ⊗ CosFourier(K, 0.5), rand((K+1) ^ 2))
# # \otimes

# newton!((F, DF, u) -> (F!(F, u, λ), DF!(DF, u, λ)), u)


# using GLMakie

# surface([x for x in LinRange(0, π, 101)], [y for y in LinRange(0, π, 101)],
#     [u(x, y) for x in LinRange(0, π, 101), y in LinRange(0, π, 101)])
