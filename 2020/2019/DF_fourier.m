function DF = DF_fourier(x,r,alpha)
m = length(x)/2;
omega = x(1);
a = x(2:end);
k = (-m+1:m-1)';
beta = (1-exp(-1i*k*omega))./(1i*k*omega);
beta(m) = 1;
beta1 = exp(-1i*k*omega)-beta;

DF = zeros(2*m);

% 1
% DF(1,1) = 0;

% 2
DF(1,2:end) = 1;

% 3
DF(2:end,1) = 1i*k.*a + (r/omega)*conv_ab(a,beta1.*a);

% 4
DF(2:end,2:end) = diag(1i*k*omega-r) + 2*r*alpha*convdiff(a,m)...
    + r*(convdiff2(a,beta,m));

end
function A = convdiff(a,m)
a_pad = [zeros(m-1,1);a;zeros(m-1,1)];
A = zeros(2*m-1);
for k=-m+1:m-1
    for j=-m+1:m-1
        A(k+m,j+m) = a_pad(k-j+2*m-1);
    end
end
end
function A = convdiff2(a,beta,m)
a_pad = [zeros(m-1,1);a;zeros(m-1,1)];
beta_pad = [zeros(m-1,1);beta;zeros(m-1,1)];
A = zeros(2*m-1);
for k=-m+1:m-1
    for j=-m+1:m-1
        A(k+m,j+m) = a_pad(k-j+2*m-1)*(beta_pad(k-j+2*m-1)+beta_pad(j+2*m-1));
    end
end
end


