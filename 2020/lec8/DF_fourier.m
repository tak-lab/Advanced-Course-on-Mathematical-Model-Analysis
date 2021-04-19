function DF = DF_fourier(x,alpha,beta,gamma,epsilon,tau)
N = (length(x)-2)/2;
omega = x(1);
a = x(2:end);
k = (-N:N).';

DF = zeros(2*N+2);

DF(1,2:end) = 1;
DF(2:end,1) = (2*omega*k.^2 + 1i*gamma*k - 1i*k*tau*beta.*exp(-1i*k*omega*tau)).*a;

[~,a2] = convp(a,2);

M = zeros(2*N+1);
% for k=-N:N
  for j=-N:N
    M(k+N+1, j+N+1) = 3*epsilon*a2(k-j+2*N+1);
  end
% end

L = diag(k.^2*omega^2+1i*gamma*k*omega + alpha + beta*exp(-1i*k*omega*tau));

DF (2:end,2:end) = L + M;