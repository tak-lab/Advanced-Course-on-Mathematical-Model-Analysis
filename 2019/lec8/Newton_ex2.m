tol = 5e-10;

x = 2*ones(3,1);x(2) = -x(2);

F = F_vec(x);

display(['Before step #1, ||F||_1 = ',num2str(norm(F,1))])
itercount = 0;

while (itercount<=100) && (norm(F,1) > tol)
    DF = DF_mat(x);
%     DF=finite_diff_DF(x);
    x = x - DF\F;
    F = F_vec(x);
    display(['After step # ',num2str(itercount+1),', ||F||_1 = ',num2str(norm(F,1))])
    itercount = itercount+1;
end

disp(x)

function F = F_vec(x)
F = zeros(3,1);
lam = x(1); x(1) = [];

F(1) = 2*x(1) - x(2) - lam*x(1);
F(2) = -x(1) + 2*x(2) - lam*x(2);
F(3) = x(1)^2 + x(2)^2 - 1;
end

function DF = DF_mat(x)
lam = x(1); x(1) = [];
DF = [ -x(1), 2-lam, -1;...
    -x(2), -1, 2-lam;...
    0, 2*x(1), 2*x(2)];
end