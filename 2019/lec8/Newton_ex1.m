tol = 1e-10;

x = 1e10;
F = Func(x);

display(['Before step #1, ||F||_1 = ',num2str(norm(F,1))])
itercount = 0;

while (itercount<=100) && (norm(F,1) > tol)
    DF = DFunc(x);
%     DF=finite_diff_DF(x);
    x = x - DF\F;
    F = Func(x);
    display(['After step # ',num2str(itercount+1),', ||F||_1 = ',num2str(norm(F,1))])
    itercount = itercount+1;
end


function y = Func(x)
y = x^2+2*x-5;
end

function dF = DFunc(x)
dF = 2*x+2;
end