%% Last Report
addpath('/Users/takitoshi/Dropbox/LaptopShare/matlab_work/mathmodel/lec7')
clear
tol = 5e-10;
alpha = 2/3*pi;
tau = 1.0;
%% Initial value for Newton method N = 61;
%
N = 61;
dde = @(t,y,z) -1*alpha*z(1)*(1+y(1));
ddehist = @(t) 0.4;
tspan = [0,400];
delay = tau;
f = dde23(dde,delay,ddehist,tspan);
%plot(f.y(1,floor(length(f.x)/2):end).',f.y(2,floor(length(f.x)/2):e nd).')
% hold on %plot(f.x(floor(length(f.x)/2):end),f.y(floor(length(f.x)/2):end)) %plot(f.x,f.y)
%
% approximate a period of solution
a = 200;
app_period = 6;
h=0.1;
f_tmp = deval(f,a+app_period/2:h:a+3*app_period/2);
find_period = abs(f_tmp - deval(f,a));
% eta0 = deval(f,a,1);
[~,ind] = min(find_period(1,:));
b = a+app_period/2 + h*(ind-1);
% transform the data into Fourier coefficient
x = ddefouriercoeffs(f,[a,b],N,1);
x0 = [2*pi/(b-a);x.'];
plot_solution(x0,3),
hold off
eta0 = real(sum(x0(2:end)));
x = x0;
%
%% Newton iteration
% tau = 2.5;
% x=x0;
F = F_fourier_kai(x,alpha,tau,eta0); %
display(['Before step #1, ||F||_1 = ',num2str(norm(F,1))])
itercount = 0;
%
while (itercount<=100) && (norm(F,1) > tol)
  DF = DF_fourier_kai(x,alpha,tau);
  % DF=finite_diff_DF(x,alpha,beta,gamma,epsilon,tau,eta0);
  x = x - DF\F;
  F = F_fourier_kai(x,alpha,tau,eta0);
  display(['After step # ',num2str(itercount+1),', ||F||_1 =',num2str(norm(F,1))])
  itercount = itercount+1;
end
%
if 2*pi/real(x(1))<0
  disp('fail')
end
L = 2*pi/real(x(1))
figure
plot_solution(x,3)

function F = F_fourier_kai(x,alpha,tau,eta0) 
N = (length(x)-2)/2;
omega = x(1);
a = x(2:end);
a3 = convp(a,3);
eta = sum(a)-eta0;
k = (-N:N)';
f = (1i*k*omega + alpha*exp(-1i*k*omega*tau).*(1+a)).*a;
F = [eta;f];
end

function DF = DF_fourier_kai(x,alpha,tau) 
N = (length(x)-2)/2;
omega = x(1);
a = x(2:end);
k = (-N:N).';
DF = zeros(2*N+2);
DF(1,2:end) = 1;
DF(2:end,1) = (1i*k-1i*k*tau*alpha.*exp(-1i*k*omega*tau).*(1+a)).*a;
L = diag(1i*k*omega + alpha*exp(-1i*k*omega*tau) + 2*alpha*exp(- 1i*k*omega*tau).*a);
DF (2:end,2:end) = L ;
end
