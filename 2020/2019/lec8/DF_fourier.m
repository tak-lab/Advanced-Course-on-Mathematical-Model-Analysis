function DF = DF_fourier(x)
m = length(x); % m = 4*N+3
N = (m-3)/4;

omega = x(1);
a = x(2:2*N+2);
b = x(2*N+3:end);
k = (-N:N)';

DF = zeros(m);

% DF(1,2:end) = ones(1,4*N+2);
DF(1,2:2*N+2) = ones(1,2*N+1);
DF(2:end,1) = [1i*k.*a; 1i*k.*b];

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

