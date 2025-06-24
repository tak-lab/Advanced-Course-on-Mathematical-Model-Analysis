m = 5 # length of Fourier
N = m - 1 # wave number
λ = 5.0 # parameter in SH
L = 6π # scale of spatial variable
ω = 2π / L # frequency of Fourier

function f(u, λ) # Vector field of SH1D
    ∂² = project(Derivative(2), space(u), space(u))
    # return λ * u - (I + ∂²)^2 * u - u^3
    return project(λ * u - (I + ∂²)^2 * u - u^3, space(u))
end

using RadiiPolynomial

function SH1D!(du, u, p, t) # Input for numerical integrator
    u_seq = Sequence(CosFourier(N, ω), u)
    du[:] = coefficients(f(u_seq, p))
end

using DifferentialEquations
u0 = randn(m)./(2.5).^(1:m)
# u0 = zeros(m)
# u0[2] = 0.2
# u0[4] = -0.15
# u0 = Sequence(CosFourier(N, ω), u0)
tspan = (0.0, 1)
p = λ
prob = ODEProblem(SH1D!, u0, tspan, p)
# sol = solve(prob)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

# using Plots
# u_seq = Sequence(CosFourier(N, ω), sol.u[150])
# plot(x -> u_seq(x),0,L/2)


using Plots
plt = plot()
for i = 1:10:length(sol.u)
    ū = Sequence(CosFourier(N, ω), sol.u[i])
    plot!(plt, x -> ū(x), 0, L, legend=false, title="Approximate solution",
        line=2,
        xlabel="\$x\$",
        ylabel="\$u\$",)
end
plot!(plt)


# Chebyshev coefficients of ū
L = 6
Ukl = zeros(m, L+1)
plt = plot()
for k = 1:m
uk(ξ) = sol( (1-ξ)*tspan[1]/2 + (1+ξ)*tspan[2]/2 )[k]
uk_seq = real(Sequence(uk, Chebyshev(L)))
Ukl[k,:] = coefficients(uk_seq)
plot!(plt,t -> abs(uk_seq(t)), -1, 1)
end
plot!(plt)


