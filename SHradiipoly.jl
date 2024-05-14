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

K = 10
λ = 1.0
u = Sequence(CosFourier(K, 0.5), rand(K+1))

newton!((F, DF, u) -> (F!(F, u, λ), DF!(DF, u, λ)), u)


using GLMakie

lines(LinRange(0, π, 101), x -> u(x); linewidth = 6)


#

_DF_ = DF!(zeros(CosFourier(K, 0.5), CosFourier(K, 0.5)), u, λ)
A = inv(_DF_)
tail_A_norm = inv(abs(λ - (1 - (K+1)^2)^2))

#

F = F!(zeros(CosFourier(3K, 0.5)), u, λ)
tail_F = copy(F)
tail_F[0:K] .= 0

Y = norm(A * F, 1) + norm(tail_F, 1) * tail_A_norm

#

DF = DF!(zeros(CosFourier(3K, 0.5), CosFourier(K, 0.5)), u, λ)

u²3 = 3u^2
tail_u²3 = copy(u²3)
tail_u²3[0:K] .= 0

Z = max(opnorm(A * DF - I, 1) + norm(tail_u²3, 1) * tail_A_norm,
        norm(u²3, 1) * inv(abs(λ - (1 - (3K+1)^2)^2)))

#

R = 2sup(Y)

W = 6max(opnorm(A, 1), tail_A_norm)*(norm(u, 1) + R)

#

interval_of_existence(interval(Y), interval(Z), interval(W), R)













# Multidim

function F!(F, u, λ)
    ∂1² = project(Derivative(2, 0), space(u), space(u))
    ∂2² = project(Derivative(0, 2), space(u), space(u))
    Δ = ∂1² + ∂2²
    project!(F, λ * u - (I + Δ)^2 * u - u^3)
    return F
end

function DF!(DF, u, λ)
    ∂1² = project(Derivative(2, 0), space(u), space(u))
    ∂2² = project(Derivative(0, 2), space(u), space(u))
    Δ = ∂1² + ∂2²
    add!(DF, λ * I - (I + Δ)^2, Multiplication(-3u^2))
    return DF
end

using RadiiPolynomial

K = 20
λ = 1.0
u = Sequence(CosFourier(K, 0.5) ⊗ CosFourier(K, 0.5), rand((K+1) ^ 2))
# \otimes

newton!((F, DF, u) -> (F!(F, u, λ), DF!(DF, u, λ)), u)


using GLMakie

surface([x for x in LinRange(0, π, 101)], [y for y in LinRange(0, π, 101)],
    [u(x, y) for x in LinRange(0, π, 101), y in LinRange(0, π, 101)])
