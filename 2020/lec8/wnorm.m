function value = wnorm(a,nu)
m = length(a); % m = 2*N+1
N = (m-1)/2;
k = (-N:N).';
w = nu.^k;
value = sum(abs(a).*w);