function value = wnorm_mat(A,nu)
m = size(A,1); % m = 2*N+1
N = (m-1)/2;
k = (-N:N).';
w = nu.^k;
value = max(sum(w.*abs(A),1)./w');
% value = sum(abs(a).*w);