function dy = func(t,y)
% N = n-1
N = length(y);
dy = zeros(N,1);

n = N+1;
dy(1) = n^2*(y(2)-2*y(1));
for k=2:N-1
  dy(k) = n^2*(y(k+1)-2*y(k)+y(k-1));
end
dy(N) = n^2*(-2*y(N)+y(N-1));
end