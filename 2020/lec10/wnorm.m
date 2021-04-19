function value = wnorm(a,nu) % the input should be column (tate) vector
% the weighted 1-norm of vector
m = length(a); % m = 2*N+1
N = (m-1)/2;
k = (-N:N).';
w = nu.^abs(k);% interval
value = sum(abs(a).*w);