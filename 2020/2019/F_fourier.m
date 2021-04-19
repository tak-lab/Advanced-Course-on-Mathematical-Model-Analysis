function F = F_fourier(x,r,alpha)
m = length(x)/2;
omega = x(1);
a = x(2:end);
k = (-m+1:m-1)';
beta = (1-exp(-1i*k*omega))./(1i*k*omega);
beta(m) = 1;

F = zeros(2*m,1);


F(1) = sum(a)-1/(1+alpha);
F(2:end) = (1i*k*omega-r).*a + r*alpha*conv_ab(a,a) + r*conv_ab(a,beta.*a);
