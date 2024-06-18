# solve Swift-Hohenberg

function F!(F, u, λ)
    ∂² = project(Derivative(2), space(u), space(u))
    project!(F, λ * u - (I + ∂²)^2 * u - u^3)
    return F
end

function DF!(DF, u, λ)
    ∂² = project(Derivative(2), space(u), space(u))
    add!(DF, λ * I - (I + ∂²)^2, Multiplication(-3u^2))
    return DF
end

using RadiiPolynomial

m = 10
N = m-1
λ = 1.0
ω = 1.0
ν = 1.0
u = Sequence(CosFourier(N, ω), rand(N+1))

newton!((F, DF, u) -> (F!(F, u, λ), DF!(DF, u, λ)), u)


# using GLMakie

# lines(LinRange(0, π, 101), x -> u(x); linewidth = 6)
using Plots
plot(x -> u(x),0 ,π, legend=false, title = "Plot cosine Fourier series",
        line=2,
        xlabel = "\$x\$",
        ylabel = "\$u(x)\$",)

# A† , An
_DF_ = DF!(zeros(CosFourier(N, ω), CosFourier(N, ω)), u, λ)
A = inv(_DF_)
tail_A_norm = inv(abs(λ - (1 - (N+1)^2)^2))

# Y0 bound
F = F!(zeros(CosFourier(3N, ω)), u, λ)
tail_F = copy(F)
tail_F[0:N] .= 0

Y = norm(A * F, Ell1(GeometricWeight(ν))) + norm(tail_F, Ell1(GeometricWeight(ν))) * tail_A_norm

#

# Z
DF = DF!(zeros(CosFourier(3N, ω), CosFourier(N, ω)), u, λ)

u²3 = 3u^2
tail_u²3 = copy(u²3)
tail_u²3[0:N] .= 0

# Z0 + Z1
Z = max(opnorm(A * DF - I, Ell1(GeometricWeight(ν))) + norm(tail_u²3, Ell1(GeometricWeight(ν))) * tail_A_norm,
        norm(u²3, Ell1(GeometricWeight(ν))) * inv(abs(λ - (1 - (3N+1)^2)^2)))

#

R = 2sup(Y)

W = 6max(opnorm(A, 1), tail_A_norm)*(norm(u, 1) + R)

#

interval_of_existence(interval(Y), interval(Z), interval(W), R)



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
