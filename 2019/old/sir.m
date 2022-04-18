function y = sir(t,x)
beta = 0.05;
nu = 0.2;
y = zeros(2,1);
y(1) = -beta*x(1)*x(2);
y(2) = beta*x(1)*x(2)-nu*x(2);