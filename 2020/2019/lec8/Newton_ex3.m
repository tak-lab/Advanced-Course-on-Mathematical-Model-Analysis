tol = 1e-12;
n = 1000;
A = randn(n);

x = randn(n+1,1);

F = F_vec(x,A);
display(['Before step #1, ||F||_1 = ',num2str(norm(F,1))])
itercount = 0;

while (itercount<=100) && (norm(F,1) > tol)
    DF = DF_mat(x,A);
%     DF=finite_diff_DF(x);
    x = x - DF\F;
    F = F_vec(x,A);
    display(['After step # ',num2str(itercount+1),', ||F||_1 = ',num2str(norm(F,1))])
    itercount = itercount+1;
end

norm(A*x(2:end)-x(1)*x(2:end),1)

function F = F_vec(x,A)
lam = x(1); x(1)=[];
F = [ x'*x-1; A*x-lam*x];
end

function DF = DF_mat(x,A)
lam = x(1); x(1)=[];
DF = [0, 2*x.';...
    -x, A-lam*eye(size(A))];
end