tol = 5e-10;
% alpha = -0.051;
% beta = -0.4;
alpha = -0.1;
beta = -5;
gamma = -2;
epsilon = -10;
tau = 2.5;
%% Initial value for Newton method
N = 61;
%
dde = @(t,y,z) [y(2); gamma*y(2) + alpha*y(1) + beta*z(1) + epsilon*y(1).^3];
ddehist = @(t) [0.4;0];
tspan = [0,400];
delay = tau;
f = dde23(dde,delay,ddehist,tspan);
plot(f.y(1,floor(length(f.x)/2):end).',f.y(2,floor(length(f.x)/2):end).')
hold on
% plot(f.y(1,:),f.y(2,:))
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
plot_solution(x0,3),hold off
eta0 = real(sum(x0(2:end)));
x = x0;
% 
%% Newton iteration
% tau = 2.5;
% x=x0;
F = F_fourier(x,alpha,beta,gamma,epsilon,tau,eta0);
% 
display(['Before step #1, ||F||_1 = ',num2str(norm(F,1))])
itercount = 0;
%
while (itercount<=100) && (norm(F,1) > tol)
    DF = DF_fourier(x,alpha,beta,gamma,epsilon,tau);
%     DF=finite_diff_DF(x,alpha,beta,gamma,epsilon,tau,eta0);
    x = x - DF\F;
    F = F_fourier(x,alpha,beta,gamma,epsilon,tau,eta0);
    display(['After step # ',num2str(itercount+1),', ||F||_1 = ',num2str(norm(F,1))])
    itercount = itercount+1;
end
% 
if  2*pi/real(x(1))<0
  disp('fail')
end
L = 2*pi/real(x(1))
figure
plot_solution(x,3)