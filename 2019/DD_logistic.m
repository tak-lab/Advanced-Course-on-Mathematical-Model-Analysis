tol = 5e-10;
m = 50;

r = 14;
alpha = 0;

%%
% N = m-1;
% dx = 1/(2*N-1);
% x = dx*(0:2*N)';
% [~,~,dn] = ellipj(x,0.25);
% x0 = fftshift(fft(dn));
% x=[5;x0];
% omega = 1;
% a = zeros(2*m-1,1);
% a(m-1) = 1;
% a(m+1) = 1;

% x = [omega;a];
F = F_fourier(x,r,alpha);

display(['Before step #1, ||F||_1 = ',num2str(norm(F,1))])
itercount = 0;

while (itercount<=100) && (norm(F,1) > tol)
    DF = DF_fourier(x,r,alpha);
%     DF=finite_diff_DF(x);
    x = x - DF\F;
    F = F_fourier(x,r,alpha);
    display(['After step # ',num2str(itercount+1),', ||F||_1 = ',num2str(norm(F,1))])
    itercount = itercount+1;
end

L = 2*pi/real(x(1))

x_norm = sum(abs(x(2:end)));
plot(r,x_norm,'.')

% plot_solution(x,1), hold on

% nu = 1.02;
% 
% [r,success] = int_radii_polynomial(u,r,alpha,nu)
