ODEformats
chebfunpref.setDefaults('factory');

%% aliasing error
f1=chebfun('sin(pi*x)');
f2=chebfun('sin(9*pi*x)');
plot(f1,LW,1.6), hold on
plot(f2,LW,1.6), ylim([-1,1])
x = -1:0.25:1; plot(x,f1(x),'r.')
hold off

%% a function f
f = @(t) exp(-(100*(t-.5).^2));
% f = @(t) exp(-erf(t.^2)+(t-1).^3);
fc = chebfun(f);
M = length(fc);
plot(fc, LW, 1.6)

%% Compute the convolution by FFT algorithm
cc = chebcoeffs(fc);% Chebyshev coefficient
a  = [flipud(0.5*cc(2:end));cc(1);0.5*cc(2:end)];% Fourier coefficient
% for p = 2:20
p = 19;
  N = (p-1)*M;
  ta = [zeros(N,1);a;zeros(N,1)];% 1. Padding zeros
  tb = ifft(ifftshift(ta));% 2. IFFT of ta
  tb2 = tb.^p;% 3. tb*^tb
  cc2 = real(fftshift(fft(tb2))*(2*p*M-1)^(p-1));% 4. FFT of tb2
  cc2_cc = [cc2(N+M);2*cc2(N+M+1:N+M+p*(M-1))];% Take Chevshev coefficient
  
  plot(chebfun(cc2_cc,'coeffs'), LW, 1.6) % If f is real function, taking real of coefficients.
%   pause
% end

%% Chebyshev coefficients by chebfun
fc2 = fc.^p;
% plot(fc2)
cfc2 = chebcoeffs(fc2);
% disp([cc2_cc(1:length(fc2)), cfc2])
norm(cc2_cc(1:length(fc2))-cfc2,1)
plotcoeffs(fc2)
hold on
plot(0:(length(cc2_cc)-1),abs(cc2_cc))
hold off

%% Application: residual of an approximate solution by Chebfun
% solve a nonlinear IVP by chebfun
tmax = 10;
p = 5;
q = 3;
N = chebop(0,tmax);
N.op = @(t,y) diff(y) - y + y.^p - y.^q;
N.lbc = 0.5;
y = N\0;
plot(y,LW,1.6)

%% calculate the DFT by FFT algorithm
yc = chebcoeffs(y);
M = length(y);
N = (p-1)*M;
ycp = convp(yc,p);
ycq = convp(yc,q);
L = length(ycp)-length(ycq);
ycq = [ycq;zeros(L,1)];

%%%%
% Differentiate the Chebyshev series by Chebfun
% dy = diff(y);
% dyc = chebcoeffs(dy);
% dyc = [dyc; zeros((p-1)*(M-1)+1,1)];
% Differentiate the Chebyshev series by coefficients
rescaleFactork = tmax/2;
dyc2 = ChebDerCoeffs(yc,0)/rescaleFactork;
dyc2 = [dyc2; zeros((p-1)*(M-1)+1,1)];
% norm(dyc-dyc2,1)
% disp([dyc,dyc2])

%% calculate the residual of ODE
yc = [yc;zeros((p-1)*(M-1),1)];
res = dyc2 - yc + ycp - ycq;
norm(res,1)

