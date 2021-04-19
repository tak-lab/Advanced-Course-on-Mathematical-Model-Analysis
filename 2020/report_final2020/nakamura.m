tol = 5e-10;
% alpha = -0.051;
% beta = -0.4;
alpha = 2;
%% Initial value for Newton method 
N = 61;
%
dde = @(t,y,z) -alpha*z(1)*(1+y(1)); ddehist = @(t) 1;
tspan = [0,10];
delay = 1;
f = dde23(dde,delay,ddehist,tspan); plot(f.x, f.y)
%%
% plot(f.y(1,:),f.y(2,:)) %
% approximate a period of solution
a = 1;
app_period = 4;
h=0.1;
f_tmp = deval(f,a+app_period/2:h:a+3*app_period/2); find_period = abs(f_tmp - deval(f,a));
% eta0 = deval(f,a,1);
[~,ind] = min(find_period(1,:));
b = a+app_period/2 + h*(ind-1);
% transform the data into Fourier coefficient 
x = ddefouriercoeffs(f,[a,b],N,1);
x0 = [2*pi/(b-a);x.'];
plot_solution(x0,3),hold off
eta0 = real(sum(x0(2:end))); x = x0;
%
%% Newton iteration
% tau = 2.5;
% x=x0;
F = F_fourier(x,alpha,eta0);
%
display(['Before step #1, ||F||_1 = ',num2str(norm(F,1))])
itercount = 0;
%
while (itercount<=100) && (norm(F,1) > tol)
  DF = DF_fourier(x,alpha);
  x = x - DF\F;
  F = F_fourier(x,alpha,eta0);
  display(['After step # ',num2str(itercount+1),', ||F||_1 = ',num2str(norm(F,1))])
  itercount = itercount+1; 
end
%
if 2*pi/real(x(1))<0
  disp('fail') 
end
L = 2*pi/real(x(1)); 
figure 
plot_solution(x,3)
% %% plot ans t = 0:0.01:20; y = 0:0.01:20;
% for ind=0:2000 
%   tmp = 0;
%   for k = -N:N
%     tmp = tmp + x(k+N+2)*exp(1i*k*x(1)*t(ind+1));
%   end
%   y(ind+1) = tmp; end
% plot(t, y)

function F = F_fourier(x,alpha,eta0) 
N = (length(x)-2)/2;
omega = x(1);
a = x(2:end);
aconv = (-N:N)'; 
for k = -N:N
  tmp = 0;
  for k_1 = -N:N
    if abs(k-k_1) > N 
      continue
    end
    tmp = tmp + exp(-1i*k_1*omega)*a(k_1+N+1)*a(k-k_1+N+1); 
  end
  aconv(k+N+1) = tmp; 
end
eta = sum(a)-eta0;
k = (-N:N)';
f = (1i*k*omega + alpha*exp(-1i*k*omega)).*a + alpha*aconv;
F = [eta;f];
end

function DF = DF_fourier(x,alpha) 
N = (length(x)-2)/2;
omega = x(1);
a = x(2:end);
aconv = (-N:N)'; 
for k = -N:N
  tmp = 0;
  for k_1 = -N:N
    if abs(k-k_1) > N 
      continue
    end
    tmp = tmp + k_1*exp(-1i*k_1*omega)*a(k_1+N+1)*a(k-k_1+N+1); 
  end
  aconv(k+N+1) = tmp; 
end
k = (-N:N).';
DF = zeros(2*N+2);
DF(1,2:end) = 1;
DF(2:end,1) = (1i*k - 1i*k*alpha.*exp(-1i*k*omega)).*a - 1i*alpha*aconv;
L = diag(1i*k*omega+alpha*exp(-1i*k*omega));
M = zeros(2*N+1); 
for k=-N:N
  for j=-N:N
    if abs(k-j) > N
      continue 
    end
    M(k+N+1, j+N+1) = alpha*(exp(-1i*j*omega) + exp(-1i*(k-j)*omega))*a(k-j+N+1); 
  end
end
DF (2:end,2:end) = L + M;
end