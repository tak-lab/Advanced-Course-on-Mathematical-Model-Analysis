function [Df] = finite_diff_DF(u,alpha,beta,gamma,epsilon,tau,eta0)

h=1e-8;
m=length(u);
E=eye(m);
Df=zeros(m);
for j=1:m
    uh = u+h*E(:,j);
    Df(:,j)=(F_fourier(uh,alpha,beta,gamma,epsilon,tau,eta0)-  F_fourier(u,alpha,beta,gamma,epsilon,tau,eta0))/h;
end
end