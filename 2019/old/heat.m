function du = heat(t,u)

n = length(u) + 1; % N = n-1 <==> n=N+1

du = zeros(n-1,1);

du(1) = n^2*(u(2)-2*u(1));
for k=2:n-2
  du(k) = n^2*(u(k+1)-2*u(k)+u(k-1));
end
du(n-1) = n^2*(-2*u(n-1)+u(n-2));

du = .01*du;