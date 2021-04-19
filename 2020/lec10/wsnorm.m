function value = wsnorm(a,nu) % the input should be the row (yoko) vector
% the norm of dual space of the weighted ell^1
m = length(a); % m = 2*N+1
N = (m-1)/2;
k = (-N:N);
w = nu.^abs(k);% interval
value = max(abs(a)./w);