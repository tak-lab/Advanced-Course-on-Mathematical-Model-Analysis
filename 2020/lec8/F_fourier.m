function F = F_fourier(x,alpha,beta,gamma,epsilon,tau,eta0)
N = (length(x)-2)/2;
omega = x(1);
a = x(2:end);
a3 = convp(a,3);
eta = sum(a)-eta0;

k = (-N:N)';
f = (k.^2*omega^2 + 1i*gamma*k*omega + alpha + beta*exp(-1i*k*omega*tau)).*a + epsilon*a3;

F = [eta;f];