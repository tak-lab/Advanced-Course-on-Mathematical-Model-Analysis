function y = sirs(t,x)
beta = 0.05;
nu = 0.2;
mu = 0.1;
y = zeros(3,1);
y(1) = -beta*x(1)*x(2) + mu*x(3);
y(2) = beta*x(1)*x(2)-nu*x(2);
y(3) = nu*x(2) - mu*x(3);