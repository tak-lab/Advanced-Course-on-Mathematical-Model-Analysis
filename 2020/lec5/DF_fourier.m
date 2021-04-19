function DF = DF_fourier(x,u0)
m = length(x); % m = 4*N+3
N = (m-3)/4;

omega = x(1);
a = x(2:2*N+2);
b = x(2*N+3:end);
k = (-N:N)';

% v0 = [3;2];
% g = @(y) [y(1)-y(1)*y(2); -y(2)+y(1)*y(2)];
% v0 = g(u0);
% bomega = bx(1);
% ba = bx(2:2*N+2);
% bb = bx(2*N+3:end);
% bomega = bx(1);
% f = f_vector_field(bx);
% af = f(2:2*N+2)/bomega;
% bf = f(2*N+3:end)/bomega;


DF = zeros(m);

% DF(1,2:end) = ones(1,4*N+2);
% DF(1,2:2*N+2) = ones(1,2*N+1);
% DF(1,2:2*N+2) = v0(1)*ones(1,2*N+1);
% DF(1,2*N+3:end) = v0(2)*ones(1,2*N+1);
% DF(1,2:2*N+2) = 1i*(-k).*flipud(ba);
% DF(1,2*N+3:end) = 1i*(-k).*flipud(bb);
% DF(1,2:2*N+2) = flipud(af);
% DF(1,2*N+3:end) = flipud(bf);
% DF(2:end,1) = [1i*k.*a; 1i*k.*b];
DF(1,1) = 0;
DF(1,2:2*N+2) = 2*(a_k_j(a,N)-u0(1));
DF(1,2*N+3:end) = 2*(a_k_j(b,N)-u0(2));


% 1
DF1 = diag(1i*k*omega-1) + convdiff(b,N);
% 2
DF2 = - convdiff(b,N);
% 3
DF3 = convdiff(a,N);
% 4
DF4 = diag(1i*k*omega+1) - convdiff(a,N);


DF(2:2*N+2,2:2*N+2) = DF1;
DF(2*N+3:end,2:2*N+2) = DF2;
DF(2:2*N+2,2*N+3:end) = DF3;
DF(2*N+3:end,2*N+3:end) = DF4;
end

function A = convdiff(b,N)
b_pad = [zeros(N,1);b;zeros(N,1)];
A = zeros(2*N+1);
for k=-N:N
    for j=-N:N
        A(k+N+1,j+N+1) = b_pad(k-j+2*N+1);
    end
end
end

function akj = a_k_j(a,N)
a_pad = [zeros(N,1);a;zeros(N,1)];
akj = zeros(2*N+1,1);
k = -N:N;
for j=-N:N
  akj(j+N+1) = sum(a_pad(k-j+2*N+1));
end
end

