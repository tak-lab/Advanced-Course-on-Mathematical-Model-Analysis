function du = advec(t,u)

N = length(u); % N = n+1 <==> n=N-1
n = N-1;

du = zeros(N,1);

x = (0:N-1)'/n*2*pi;
cx = coef(x);

% du(1) = -n*(u(1)-u(N-1))*cx(1);
du(1) = -n*(u(2)-u(N-1))*cx(1)/2;
for k=2:N-1
%   du(k) = -n*(u(k)-u(k-1))*cx(k);
  du(k) = -n*(u(k+1)-u(k-1))*cx(k)/2;
end
% du(N) = -n*(u(N)-u(N-1))*cx(N);
du(N) = -n*(u(2)-u(N-1))*cx(N)/2;
end

function c = coef(x)
% c = 0.2*ones(size(x));
c = 0.2 + sin(x-1).^2;
end