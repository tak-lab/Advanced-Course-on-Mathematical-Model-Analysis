% MATLABÇÃê⁄ë±ï˚ñ@
% https://www.u.tsukuba.ac.jp/remote/
% https://www.u.tsukuba.ac.jp/remote/#view2win

%% SIR
[t,y] = ode45(@(t,x) sir(t,x),[0,5],[10,1]);
plot(t,y)
hold on
plot(t,11-y(:,1)-y(:,2))

%% SIRS
[t,y] = ode45(@(t,x) sirs(t,x),[0,60],[10,1,0]);
plot(t,y)

%% Lotka-Volterra
a=1;b=1;c=0.2;d=1;
[t,y] = ode45(@(t,x) LotkaVolterra(t,x,a,b,c,d),[0,60],[1,1]);
plot(t,y)
plot(y(:,1),y(:,2))
%
y = chebfun.ode45(@(t,x) LotkaVolterra(t,x,a,b,c,d),[0,60],[1,1]);
plot(y(:,1),y(:,2))
%
tol=1e-20;opts = odeset('abstol',tol,'reltol',tol);
y = chebfun.ode45(@(t,x) LotkaVolterra(t,x,a,b,c,d),[0,60],[1,1],opts);
plot(y(:,1),y(:,2))

%% Chebyshev points
LW = 'linewidth'; FS = 'fontsize'; MS = 'markersize'; ms = 14; lw = 0.7;
IN = 'interpret'; LT = 'latex';
n = 40;
j = (0:n)';
thetaj = j*pi/n;
zj = exp(1i*thetaj);
xj = real(zj);
hold off, plot(zj,'.-k'),ylim([0 1.1]), xlim([-1,1])
title('Equispaced points on the unit circle')

hold on
for j = 2:n
  plot([xj(j) zj(j)],'k','linewidth',0.7)
end
plot(xj,0*xj,'.r'), title('Chebyshev points')
hold off

%% Chebyshev polynomial
f=chebpoly(n);
hold on
plot(f)
hold off

%% Function values on chebyshev points using chebfun
% domain = [0,1];%r=.5;cent=.5;
% fc = chebfun('exp(erf(x^2)+x.^5).*sin(5*pi*x) + x',domain);
% n = length(fc) - 1;
% plot(fc), hold on
% cpts = chebpts(n+1,domain);
% fvals=fc(cpts);
% plot(cpts,fvals,'.',MS,14)
% title('$f(x)=\exp(\mathrm{erf}(x^2)+x^5)\sin(5\pi x)+x$',IN,LT)
% grid on, hold off, xlim(domain)

%% Function values on chebyshev points without using chebfun
a = 0; b = 1;
domain = [a,b];
fc = @(x) exp(erf(x.^2)+x.^5).*sin(5*pi*x) + x;
fplot(fc,domain)
hold on
xj = 0.5*(1-xj)*a+0.5*(1+xj)*b;
fvals=fc(xj);
plot(xj,fvals,'.',MS,14)
title('$f(x)=\exp(\mathrm{erf}(x^2)+x^5)\sin(5\pi x)+x$',IN,LT)
grid on, hold off, xlim(domain)

%% Convert function values into Chebyshev coefficients
valsUnitDisc = [flipud(fvals); fvals(2:end-1)];
FourierCoeffs = real(fft(valsUnitDisc));
ChebCoeffs = (FourierCoeffs(1:n+1))/n;
ChebCoeffs(1) = ChebCoeffs(1)/2;
ChebCoeffs(end) = ChebCoeffs(end)/2;
% ChebCoeffs2 = chebvals2chebcoeffs(fvals);
% ChebCoeffs3 = chebcoeffs(fc);
% g = real(fft(fvals([1:n+1 n:-1:2])))/(2*n); % FFT
% a = [g(1); g(2:n)+g(2*n:-1:n+2); g(n+1)];
% format long
semilogy(abs(ChebCoeffs))
% display([ChebCoeffs,ChebCoeffs2,ChebCoeffs3])
% display([chebpoly(fc)' ChebCoeffs chebpoly(fc)'-ChebCoeffs])

