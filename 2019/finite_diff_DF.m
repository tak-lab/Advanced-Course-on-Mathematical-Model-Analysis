function [Df] = finite_diff_DF(u)

h=1e-8;
m=length(u);
E=eye(m);
Df=zeros(m);
for j=1:m
    uh = u+h*E(:,j);
    Df(:,j)=( F_fourier(uh)- F_fourier(u))/h;
end
end