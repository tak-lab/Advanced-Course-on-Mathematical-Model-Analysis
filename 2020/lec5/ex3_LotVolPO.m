clf
clear
tol = 5e-10;
N = 71;
u0 = [4;1];

%% Initial guess by ode45
tspan = [0,100];
opts = odeset('RelTol',1e-14);
f = ode45(@(t,y) [y(1)-y(1)*y(2); -y(2)+y(1)*y(2)],tspan,u0,opts);
plot(f.y(1,floor(length(f.x)/2):end).',f.y(2,floor(length(f.x)/2):end).'),hold on
% % plot(f.y(1,:),f.y(2,:))
% approximate a period of solution
a = 0;
app_period = 8;
h=0.1;
f_tmp = deval(f,a+app_period/2:h:a+3*app_period/2);
find_period = abs(f_tmp - deval(f,a));
eta0 = deval(f,a,1);
[~,ind] = min(find_period(1,:));
b = a+app_period/2 + h*(ind-1);
% transform the data into Fourier coefficient
a_fourier = fouriercoeffs(f,[a,b],N,1);
b_fourier = fouriercoeffs(f,[a,b],N,2);

%% Newton iteration
omega = 2*pi/(b-a);
x0 = [omega;a_fourier.';b_fourier.'];
x = x0;

% x0 = randn(4*N+3,1);
% k = (-N:N)';
% s = 1.5;
% a0 = randn./k.^s; a0(N+1) = randn;
% b0 = randn./k.^s; b0(N+1) = randn;
% x0 = [randn; a0; b0];
% x = x0;


F = F_fourier(x,u0);

display(['Before step #1, ||F||_1 = ',num2str(norm(F,1))])
itercount = 0;

while (itercount<=10) && (norm(F,1) > tol)
    DF = DF_fourier(x,u0);
%     DF=finite_diff_DF(x);
    x = x - DF\F;
    F = F_fourier(x,u0);
    display(['After step # ',num2str(itercount+1),', ||F||_1 = ',num2str(norm(F,1))])
    itercount = itercount+1;
end

L = 2*pi/real(x(1))
plot_solution(x,3)
DF = DF_fourier(x,u0);
cond(DF)

% nu = 1.02;
% 
% [r,success] = int_radii_polynomial(u,r,alpha,nu)
