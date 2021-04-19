%% 
ODEformats
%% Chebyshev points on the unit circle
n = 10;
tt = linspace(0,pi,n+1);
zz = exp(1i*tt);
hold off, plot(zz,'.-k'), axis equal, ylim([0 1.1])%, xlim([-1,1])
title('Equispaced points on the unit circle')

%% Chebyshev points on the real line
xx = chebpts(n+1);
hold on
for j = 2:n
plot([xx(n+2-j) zz(j)],'k',LW,0.7)
end
plot(xx,0*xx,'.r'), title('Chebyshev points')
hold off

%% Chebyshev polynomials of degree n
f=chebpoly(n);
hold on
plot(f)
hold off


%% Function values on the Chebyshev points
domain = [0,1];r=.5;cent=.5;
fc = chebfun('exp(erf(x^2)+x.^5).*sin(5*pi*x) + x',domain);
n = length(fc) - 1;
% LW = 'linewidth'; FS = 'fontsize'; MS = 'markersize'; ms = 14; lw = 0.7;
plot(fc), hold on
cpts = chebpts(n+1,domain);
fvals=fc(cpts);
plot(cpts,fvals,'.',MS,14)
title('$f(x)=\exp(\mathrm{erf}(x^2)+x^5)\sin(5\pi x)+x$',IN,LT)
grid on, hold off, xlim(domain)


%% Convert the values on the Chebyshev points to Chebyshev coefficients
valsUnitDisc = [flipud(fvals); fvals(2:end-1)];
FourierCoeffs = real(fft(valsUnitDisc));
ChebCoeffs = (FourierCoeffs(1:n+1))/n;
ChebCoeffs(1) = ChebCoeffs(1)/2;
ChebCoeffs(end) = ChebCoeffs(end)/2;
ChebCoeffs2 = chebvals2chebcoeffs(fvals);
ChebCoeffs3 = chebcoeffs(fc);
% g = real(fft(fvals([1:n+1 n:-1:2])))/(2*n); % FFT
% a = [g(1); g(2:n)+g(2*n:-1:n+2); g(n+1)];
format long
display([ChebCoeffs,ChebCoeffs2,ChebCoeffs3])
% display([chebpoly(fc)' ChebCoeffs chebpoly(fc)'-ChebCoeffs])