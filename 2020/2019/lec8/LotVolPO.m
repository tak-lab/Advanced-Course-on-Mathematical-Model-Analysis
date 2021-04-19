tol = 5e-10;
N = 100;
% %%
L = chebop(0,10); L.lbc = [4;1];
L.op = @(t,u,v) [diff(u)-u+u*v; diff(v)+v-u*v];
[x,y]=L\0;%plot(x,y)
tx = chebcoeffs(x);
a = [flipud(.5*tx(2:N));tx(1);.5*tx(2:N)];
ty = chebcoeffs(y);
b = [flipud(.5*ty(2:N));ty(1);.5*ty(2:N)];
% a=a/2;
% b=b/2;

%%
omega = 10;
% a = 1e-5*(rand(2*N+1,1)+rand(2*N+1,1)*1i);
% % a(N+1) = .5;
% a(N:N+2) = .5*[-1,2,-1];
% b = 1e-5*(rand(2*N+1,1)+rand(2*N+1,1)*1i);
% b(N:N+2) = 5*[-1,2,-1];

x = [omega;a;b];
F = F_fourier(x);

display(['Before step #1, ||F||_1 = ',num2str(norm(F,1))])
itercount = 0;

while (itercount<=100) && (norm(F,1) > tol)
    DF = DF_fourier(x);
%     DF=finite_diff_DF(x);
    x = x - DF\F;
    F = F_fourier(x);
    display(['After step # ',num2str(itercount+1),', ||F||_1 = ',num2str(norm(F,1))])
    itercount = itercount+1;
end

L = 2*pi/real(x(1))
% plot_solution(u,2)

% nu = 1.02;
% 
% [r,success] = int_radii_polynomial(u,r,alpha,nu)
