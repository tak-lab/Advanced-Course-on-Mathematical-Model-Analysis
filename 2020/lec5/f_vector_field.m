function F = f_vector_field(x)
m = length(x); % m = 4*N+3
N = (m-3)/4;
F = zeros(m,1);

omega = x(1);
a = x(2:2*N+2);
b = x(2*N+3:end);
k = (-N:N)';


F(1) = 0;

F(2:2*N+2) = a - conv_ab(a,b);
F(2*N+3:end) = - b + conv_ab(a,b);
% F(1) = sum(a)+sum(b)-5;

