function du = fujita_neu(t,u)

N = size(u,1);
du = zeros(N,1);
n = N+1;
p = 3; % exponent of reaction term

du(1) = n^2*(u(2)-u(1))+u(1)^p;
for k=2:n-2
  du(k) = n^2*(u(k+1)-2*u(k)+u(k-1))+u(k)^p;
end
du(n-1) = n^2*(-u(n-1)+u(n-2))+u(n-1)^p;