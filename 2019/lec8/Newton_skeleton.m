tol = 5e-10;
F = F(x);

display(['Before step #1, ||F||_1 = ',num2str(norm(F,1))])
itercount = 0;

while (itercount<=100) && (norm(F,1) > tol)
    DF = DF_fourier(x);
%     DF=finite_diff_DF(x);
    x = x - DF\F;
    F = F(x);
    display(['After step # ',num2str(itercount+1),', ||F||_1 = ',num2str(norm(F,1))])
    itercount = itercount+1;
end
