function F = F_fourier(x,u0)
m = length(x); % m = 4*N+3
N = (m-3)/4;
F = zeros(m,1);

omega = x(1);
a = x(2:2*N+2);
b = x(2*N+3:end);
k = (-N:N)';

% v0 = [3;2];
% u0 = [4;1];
% g = @(y) [y(1)-y(1)*y(2); -y(2)+y(1)*y(2)];
% v0 = g(u0);
% F(1) = -u0'*v0 + (v0(1)*sum(a) + v0(2)*sum(b));

% bomega = bx(1);
% f = f_vector_field(bx);
% af = f(2:2*N+2)/bomega;
% bf = f(2*N+3:end)/bomega;

% F(1) = sum(flipud(a).*(1i.*k.*ba)) + sum(flipud(b).*(1i.*k.*bb));
% F(1) = sum(flipud(a).*(af)) + sum(flipud(b).*(bf));
F(1) = (sum(a)-u0(1))^2 + (sum(b)-u0(2))^2;

% F(1) = sum(a)-4;
F(2:2*N+2) = 1i*k*omega.*a - a + conv_ab(a,b);
F(2*N+3:end) = 1i*k*omega.*b + b - conv_ab(a,b);
% F(1) = sum(a)+sum(b)-5;

