function du = gelfand_per(t,u)

n = length(u); % N = n-1 <==> n=N+1
N = n-1;
p = 3; % the exponent of the reaction term

du = zeros(n,1);

du(1) = N^2*(u(2)-2*u(1)+u(n-1)) + exp(u(1));
for k=2:n-1
  du(k) = N^2*(u(k+1)-2*u(k)+u(k-1)) + exp(u(k));
end
du(n) = N^2*(u(2)-2*u(n)+u(n-1)) + exp(u(n));
