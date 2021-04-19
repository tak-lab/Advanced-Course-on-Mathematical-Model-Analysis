clf
clear
tol = 5e-10;
N = 71;
etaeta = 1.2;

addpath('/Users/takitoshi/Dropbox/LaptopShare/matlab_work/mathmodel/lec5')
addpath('/Users/takitoshi/Dropbox/LaptopShare/matlab_work/mathmodel/lec7')


%% Initial guess by ode45
tspan = [0,400];
opts = odeset('RelTol',1e-14);
f = ode45(@(t,y) [y(2);-y(1)^3+0.5*y(1)^2],tspan,[etaeta,0],opts);
plot(f.y(1,floor(length(f.x)/2):end).',f.y(2,floor(length(f.x)/2):end).')
%plot(f.x(1,:),f.y(1,:))
% approximate a period of solution
a = 24;
app_period = 8;
h=0.1;
f_tmp = deval(f,a+app_period/2:h:a+3*app_period/2);
find_period = abs(f_tmp - deval(f,a));
%eta0 = deval(f,a,1);
[~,ind] = min(find_period(1,:));
b = a+app_period/2 + h*(ind-1);
% transform the data into Fourier coefficient
a_fourier = fouriercoeffs(f,[a,b],N,1);
eta0 = real(sum(a_fourier));

%% Newton iteration
omega = 2*pi/(b-a);
x = [omega; a_fourier.'];

% x0 = randn(4*N+3,1);
% k = (-N:N)';
% s = 1.5;
% a0 = randn./k.^s; a0(N+1) = randn;
% b0 = randn./k.^s; b0(N+1) = randn;
% x0 = [randn; a0; b0];
% x = x0;

F = F_fourier(x,eta0);

display(['Before step #1, ||F||_1 = ',num2str(norm(F,1))])
itercount = 0;

while (itercount<=100) && (norm(F,1) > tol)
  DF = DF_fourier(x);
  %     DF=finite_diff_DF(x);
  x = x - DF\F;
  F = F_fourier(x,eta0);
  display(['After step # ',num2str(itercount+1),', ||F||_1 = ',num2str(norm(F,1))])
  itercount = itercount+1;
end

L = 2*pi/real(x(1))
plot_solution(x,2)
%plot_solution(x,3)


% nu = 1.02;
%
% [r,success] = int_radii_polynomial(u,r,alpha,nu)


function F = F_fourier(x,eta_0)
m = length(x); % m = 2*N+2
N = (m-2)/2;
F = zeros(m,1);

omega = x(1);
a = x(2:end);
%a5 = convp(a,5);
a3 = convp(a,3);
a2 = convp(a,2);
eta = sum(a)-eta_0;
k = (-N:N)';

F(1) = eta;
%F(2:end) = (-k.^2*omega^2).*a +a5 -4*a3;
F(2:end) = (-k.^2*omega^2).*a +a3 -0.5*a2;
%F(2:end) = (-k.^2*omega^2).*a +a5;
% F(1) = sum(a)+sum(b)-5;

end

function DF = DF_fourier(x)
m = length(x); % m = 2*N+2
N = (m-2)/2;

omega = x(1);
a = x(2:end);
k = (-N:N)';

DF = zeros(m);

DF(1,2:end) = ones(1,2*N+1);
DF(2:end,1) = (-2*k.^2*omega).*a;

%[~,a4] = convp(a,4);
[~,a2] = convp(a,2);

M = zeros(2*N+1);
for kk=-N:N
  for j=-N:N
    %M(k+N+1,j+N+1) = 5*a4(k-j+2*N+1) -3*4*a2(k-j+2*N+1);
    %M(kk+N+1,j+N+1) = 5*a4(kk-j+2*N+1);
    M(kk+N+1,j+N+1) = 3*a2(kk-j+2*N+1);
  end
end
M = M -2*0.5*convdiff(a,N);

L = diag(-k.^2*omega^2);

DF(2:end,2:end) = M+L;
end

function A = convdiff(b,N)
b_pad = [zeros(N,1);b;zeros(N,1)];
A = zeros(2*N+1);
for k=-N:N
  for j=-N:N
    A(k+N+1,j+N+1) = b_pad(k-j+2*N+1);
  end
end
end
