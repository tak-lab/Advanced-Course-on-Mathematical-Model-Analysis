function value = wnorm_mat(A,nu)% the weighted 1-norm of matrix
m = size(A,1); % m = 2*N+1
N = (m-1)/2;
k = (-N:N)';
w = nu.^abs(k); % interval
value = max(sum(w.*abs(A),1)./w');
