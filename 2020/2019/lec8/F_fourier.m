function F = F_fourier(x)
m = length(x); % m = 4*N+3
N = (m-3)/4;
F = zeros(m,1);

omega = x(1);
a = x(2:2*N+2);
b = x(2*N+3:end);
k = (-N:N)';

F(1) = sum(a)-0.25;
F(2:2*N+2) = 1i*k*omega.*a - a + conv_ab(a,b);
F(2*N+3:end) = 1i*k*omega.*b + b - conv_ab(a,b);
% F(1) = sum(a)+sum(b)-1/2;
