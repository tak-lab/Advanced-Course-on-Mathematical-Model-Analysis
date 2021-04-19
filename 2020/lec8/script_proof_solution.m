tol = 5e-10;
% alpha = -0.051;
% beta = -0.4;
alpha = -0.1;
beta = -5;
gamma = -2;
epsilon = -10;
tau = 2.5;
%% Initial value
N = 61;
%
dde = @(t,y,z) [y(2); gamma*y(2) + alpha*y(1) + beta*z(1) + epsilon*y(1).^3];
ddehist = @(t) [0.4;0];
tspan = [0,400];
delay = tau;
f = dde23(dde,delay,ddehist,tspan);
% plot(f.y(1,floor(length(f.x)/2):end).',f.y(2,floor(length(f.x)/2):end).'), hold on
% plot(f.y(1,:),f.y(2,:))

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
% plot_solution(x0,3),hold off
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
% figure
plot_solution(x,3)


%% verify the solution
nu = 1.1;
DF = DF_fourier(x,alpha,beta,gamma,epsilon,tau);
A = inv(DF);
% Y0
delta = F_fourier_ext(x,alpha,beta,gamma,epsilon,tau,eta0);
delta0 = delta(1);
delta1 = delta(2:end);
delta1_N = delta1(2*N+1:end-2*N);
delta1(2*N+1:end-2*N) = [];
delta1_tail = delta1;

mu_k = @(k,omega) k.^2*omega^2 + 1i*k*omega + alpha + beta*exp(-1i*k*omega*tau);

k_tail = ([-3*N:-N-1,N+1:3*N])';

Y0 = max(abs(A(1,1)*delta0 + A(1,2:end)*delta1_N),...
  wnorm(A(2:end,1)*delta0,nu) + wnorm(A(2:end,2:end)*delta1_N,nu) + wnorm(delta1_tail./abs(mu_k(k_tail,x(1))),nu))

%
% Z0
k = (-N:N)';
w = nu.^k;
B = eye(size(A)) - A*DF;
Z0_0  = abs(B(1,1)) + max(abs(B(1,2:end)).*w');
Z0_1 = wnorm(B(2:end,1),nu) + wnorm_mat(B(2:end,2:end),nu);

Z0 = max(Z0_0,Z0_1)
%% Z1
% [~,a2] = convp(x(2:end),2);
% zeta = zeros(2*N+1,1);
% for k=-N:N
%   j = k-N:-N-1; wj = nu.^j;
%   zeta1 = max(abs(a2(k-j+2*N+1))./wj);
%   j = N+1:k+N; wj = nu.^j;
%   zeta(k) = max(abs(a2(k-j+2*N+1))./wj)
% end
% wn = nu^(N+1);
% abs(B(1,1))/wn + abs(B(1,2:end))*zeta
Z1 = 0;

%% Z2
lam = x(1); x(1) = [];
% B = @(r) [0,0,0,0;...
%   (1+(2*x(3)+r)), 0,0,(1+(x(3)+r))+(lam+r);...
%   (1+(2*x(1)+r)),(1+(x(1)+r))+(lam+r),0,0;...
%   (1+(2*x(2)+r)),0, (1+(x(2)+r))+(lam+r),0];
Z2_3 = 10;
Z2_2 = 10;
Z2 = @(r) Z2_3*r^2 + Z2_2*r;
% 
p = @(r) Z2(r)*r - (1-Z1-Z0)*r + Y0;
r = fzero(p,0.01)
% disp(r)
figure
fplot(p,[0,0.12],'linewidth',1.6), hold on,
plot([0,0.12], [0,0], 'k-.')
title('The radii polynomial')
xlabel('$r$','interpreter','latex'), ylabel('$p(r)$','interpreter', 'latex')

