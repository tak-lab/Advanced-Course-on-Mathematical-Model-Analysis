using FFTW

# convolution
function quadratic(a1,a2)
    m=Int((length(a1)+1)/2)
    ta1=[zeros(m);a1;zeros(m)]; tu1=ifft(ifftshift(ta1))
    ta2=[zeros(m);a2;zeros(m)]; tu2=ifft(ifftshift(ta2))
    F=fftshift(fft(tu1.*tu2))
    return (4*m-1)*F[m+1:3*m-1]
end

# verifyfft
using IntervalArithmetic
function verifyfft(z::Vector{T}, sign=1) where T
    n = length(z); col = 1; array1 = true
    if n==1
        Z = map(T,z)
        return Z
    else
        isrow_ = false
    end
    log2n = Int(round(log2(n))) #check dimension
    if 2^log2n ≠ n #2の倍数でない場合はエラー表示
        error("length must be power of 2")
    end
    #bit-reversal(ビットリバース)
    f = 2^(log2n-1)
    v = [0;f]
    for k = 1:log2n-1
#         f = 0.5*f
        f = f >> 1
        v = append!(v,f.+v)
    end
    z2 = zeros(n,col)
    if isa(real(z[1]),Interval)
        z2 = map(T,z2)
    end
    #zを入れ替え
    for j = 1: n
        z2[j,:] = z[v[j]+1,:]
    end
    #Danielson-Lanczos algorithm
    Z = complex(map(Interval,z2))
    Index = reshape([1:n*col;],n,col)

    theta = map(Interval,sign * (0:n-1)/n); # division exact because n is power of 2
    Phi = cospi.(theta) + im*sinpi.(theta) # SLOW?

    v = [1:2:n;]
    w = [2:2:n;]
    t = Z[w,:]
    Z[w,:]  = Z[v,:] - t
    Z[v,:]  = Z[v,:] + t
    for index　in 1: (log2n-1)    
        m = 2^index
        m2 = 2*m
        vw = reshape([1:n;],m2,Int(n/m2))
        v = vw[1: m, :]
        w = vw[m+1: m2, : ]
        indexv = reshape(Index[v[:],:],m,Int(col*n/m2))
        indexw = reshape(Index[w[:],:],m,Int(col*n/m2))
        Phi1 = repeat(Phi[1:Int(n/m):end],outer=[1,Int(col*n/m2)])
        t = Phi1 .*  Z[indexw]
        Z[indexw] = Z[indexv] - t 
        Z[indexv] = Z[indexv] + t
    end
    reverse(Z[2:end,:],dims=2)
     if sign==-1
        Z = Z/n
    end
    if isrow_
        Z = transpose(Z)　#転置
    end
    if array1
        Z = Z[:,1]
    end
    return Z
end

function int_quadratic_sumFFT(a1,a2) # the convolution term using FFT
    M = length(a1)
    a1 = [reverse(a1[2:end]);a1]
    a2 = [reverse(a2[2:end]);a2]
    i2=interval(2); i4=interval(4)

    # We make sure that the inputs are powers of 2 %
    M1 = Int((2^ceil(log(4*M-1)/log(2))-4*M)/2)   # Hence 4*M+2*M1 is a power of 2 %

    ta1=[zeros(Complex{Interval{Float64}}, M+M1+1);a1;zeros(Complex{Interval{Float64}}, M+M1)];
    ta2=[zeros(Complex{Interval{Float64}}, M+M1+1);a2;zeros(Complex{Interval{Float64}}, M+M1)];

    tu1=verifyfft(ifftshift(ta1),-1)
    tu2=verifyfft(ifftshift(ta2),-1)
    F=fftshift(verifyfft(tu1.*tu2,1))
    return (i4*M+i2*M1)*F[M1+2*M+1:M1+3*M]

end


# Vector field
function F_steady_states(a) # Input: cosine Fourier
    N = length(a)-1 # maximum wave #
    omega = 2*pi
    k = (0:N)
    a_ext = [reverse(a[2:end]); a]
    a2 = quadratic(a_ext,a_ext); 
    return -k.^2*omega^2 .* a + a2[N+1:end]
end

function iF_steady_states(a) # Input: cosine Fourier
    N = length(a)-1 # maximum wave #
    omega = 2*pi
    k = (0:N)
    # a_ext = [reverse(a[2:end]); a]
    # a2 = quadratic(a_ext,a_ext);
    a2 = int_quadratic_sumFFT(a,a);
    return -k.^2*omega^2 .* a + a2
end

# Jacobian matrix
function DF_steady_states(a)
    N = length(a)-1
    omega = 2*pi
    a = [reverse(a[2:end]);a;zeros(N)]
    DF_nonlinear = zeros(ComplexF64,N+1,N+1)
    DF_nonlinear[:,1] = 2*a[N+1:2*N+1]
    n = (0:N)
    for ell = 1:N
        DF_nonlinear[:,ell+1] = 2*(a[n.+(-ell+N+1)]+a[n.+(ell+N+1)])
    end 
    return diagm(-(n.^2)*omega^2) + DF_nonlinear
end

function iDF_steady_states(a)
### INPUT for F_eigenpairs: x = (lambda,a,b) in C^{2m+1}
### a=(a0,a1,a2,...,a_{m-1}): vector of Fourier coefficients of the steady state
    m = length(a)
    k = (0:m-1)
    omega = 2*@interval(pi)
    mu = interval.(-k.^2) * omega^2
    a = map(Interval,[reverse(a[2:end]);a;zeros(ComplexF64,m-1)])
    Dfa_nonlinear = zeros(Complex{Interval{Float64}},m,m)
    Dfa_nonlinear[:,1] = 2*a[m:2*m-1];
    for ell = 1:m-1
        Dfa_nonlinear[:,ell+1] = 2*(a[k.+(-ell.+m)]+a[k.+(ell+m)])
    end 
    return diagm(mu) + Dfa_nonlinear
end


using LinearAlgebra

function norma(a,nu)
    m = length(a)
    nu_power = [1; 2*nu.^interval.(abs.(1:m-1))]
    return sum(abs.(a).*nu_power)
end

function mnorma(N,M,nu)
    nu_power = [1; 2*nu.^interval.(abs.(1:N))]
    return interval(maximum(mag,sum(abs.(M) .* nu_power, dims=1)./nu_power'));
end

function operator_norm(C,inu,tail)
    m = size(C,1)
    omega = 2 * @interval(pi);
    mu_m = -m.^2 * omega^2;
    return max(sup(mnorma(m-1,C,inu)), sup(isequal(tail,1)/abs(mu_m)))
end


# Newton-iteration
function newton(a)
    tol = 5e-11; # tolerance for Newton's method
    F = F_steady_states(a)
    nF = norm(F,1)
    println("Before iteration: $(nF)")
    k=0
    while (k<=60) && (nF > tol)
        DF = DF_steady_states(a)
        a = a - DF\F
        F = F_steady_states(a)
        nF = norm(F)
        println("After $(k+1) th iteration: $(nF)")
        k = k+1
    end
    return a
end

# Plot
using Plots
function plot_periodic_complex(a)
    N = (length(a)-1)
    omega = 2*pi
    x = (0:.001:1)
    n = length(x)
    k = (-N:N)
    a = [reverse(a[2:end]); a];
    u = zeros(ComplexF64,n)    
    for j = 1:n
        u[j] = sum(a.*exp.(im*k*omega*x[j]))
    end
    plot(x,real(u),line=2,xlabel = "\$x\$",ylabel = "\$Re(u)\$, \$Im(u)\$",label="Re(u)",)
    plot!(x,imag(u),line=2,label="Im(u)")
end


# radii-polynomial
# using IntervalArithmetic
function  rad_poly_eigenpairs(a,nu)
    success = 0
    DF = DF_steady_states(a)
    A = inv(DF)
    iA = map(Interval,A)
    m = length(a)

    ### Y0 ###
    a_ext = [a;zeros(ComplexF64,m-1)] 
    inu = interval(nu)
    F1_ext = iF_steady_states(a_ext)
    omega = 2*@interval(pi)
    k_tail = (m:2*m-2)
    mu_k_tail = omega^2 .* interval.(-k_tail.^2)
    y1_F = abs.(iA*F1_ext[1:m])#y_F(2:m+1); 
    y1_tail = abs.(F1_ext[m+1:2*m-1]./mu_k_tail)
    y1 = [y1_F;y1_tail]
    # nu_power = [1; 2*inu.^interval.(1:2*m-2)]
    # Y0 = mag(sum(y1.*nu_power))
    Y0 = norma(y1,nu)

    ### Z0 ###
    iDF = iDF_steady_states(a)
    B = I - iA*iDF
    Z0 = operator_norm(B,inu,0)

    ### Z1 ###
    hQa = zeros(m,1)#; hQb = zeros(m,1);
    for k=1:m-1
        k_prime = (m:k+m)
        hQa[k+1]=(1/2)*maximum(mag, abs.( a_ext[abs.(k.-k_prime) .+ 1] )./ (inu.^interval.(abs.(k_prime))))
    end

    hQa = interval.(hQa)
    Z1 = mag(2*norma(iA*hQa,inu) + norma(a,inu)/(2*m^2*@interval(pi)^2))

    ### Z2 ###
    normA = operator_norm(iA,inu,1)
    Z2 = 2*normA

    println("Y0 = $Y0, Z0 = $Z0, Z1 = $Z1, Z2 = $Z2")

    if inf(1-Z0-Z1)>0 
      if inf((1-Z0-Z1)^2-4*Y0*Z2) > 0  
        rmin = sup(((1-Z0-Z1) - sqrt((1-Z0-Z1)^2-4*Y0*Z2))/(2*Z2))
        rmax = inf(((1-Z0-Z1) + sqrt((1-Z0-Z1)^2-4*Y0*Z2))/(2*Z2))
        if rmin<rmax 
          success=1
          r0 = [rmin rmax]
        else
          println("failure: rmin > rmax")
        end
      else
        println("failure: discriminant is negative")
      end
    else
        println("failure: linear term is positive")
    end
    return r0, success
end