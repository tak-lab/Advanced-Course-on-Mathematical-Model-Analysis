function F = F_fourier_ext(x,alpha,beta,gamma,epsilon,tau,eta0)
N = (length(x)-2)/2;
omega = x(1);
a = [zeros(2*N,1);x(2:end);zeros(2*N,1)];
[~,a3] = convp(x(2:end),3);
eta = sum(a)-eta0;

k = (-3*N:3*N)';
f = (k.^2*omega^2 + 1i*gamma*k*omega + alpha + beta*exp(-1i*k*omega*tau)).*a + epsilon*a3;

F = [eta;f];