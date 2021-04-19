clear
tol = 5e-10;
% alpha = -0.051;
% beta = -0.4;
alpha = -0.1;
beta = -5;
gamma = -2;
epsilon = -10;
tau = 2.5;
%% Initial value
N = 53; % should be odd number
% N = 50;
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
plot_solution(x,2)


%% verify the solution
x(3:2:end)=0; % symmetry condition: a_{2*m+1} = 0
bomega = real(x(1));
ba = x(2:end);
nu = 1.1;
DF = DF_fourier(x,alpha,beta,gamma,epsilon,tau);% interval
A = inv(DF);
A_a0 = A(1,2:end);
A_a1 = A(2:end,2:end);
A_o1 = A(2:end,1);
% Y0
delta = F_fourier_ext(x,alpha,beta,gamma,epsilon,tau,eta0);% interval
delta0 = delta(1);
delta1 = delta(2:end);
delta1_N = delta1(2*N+1:end-2*N);
delta1(2*N+1:end-2*N) = 0;
delta1_tail = delta1;

mu_k = @(k,omega) k.^2*omega^2 + 1i*k*omega + alpha + beta*exp(-1i*k*omega*tau);

% k_tail = ([-3*N:-N-1,N+1:3*N])';
k_tail = (-3*N:3*N)';

Y0 = max(abs(A(1,1)*delta0 + A_a0*delta1_N),...
  wnorm(A_o1*delta0 + A_a1*delta1_N,nu)...
  + wnorm(delta1_tail./mu_k(k_tail,bomega),nu))

%
% Z0
% k = (-N:N)';
% w = nu.^abs(k);
B = eye(size(A)) - A*DF;
% Z0_0  = abs(B(1,1)) + max(abs(B(1,2:end))./w');
Z0_0  = abs(B(1,1)) + wsnorm(B(1,2:end),nu);
Z0_1 = wnorm(B(2:end,1),nu) + wnorm_mat(B(2:end,2:end),nu);

Z0 = max(Z0_0,Z0_1)

% Z1
[~,a2] = convp(ba,2);
zeta = zeros(2*N+1,1);
for k=-N:N
  j = (k-2*N:-N-1)';
  if isempty(j)
    zeta1 = -1;
  else
    wj = nu.^abs(j);
    zeta1 = max(abs(a2(k-j+2*N+1))./wj);
  end
  j = (N+1:k+N)';
  if isempty(j)
    zeta2 = -1;
  else
    wj = nu.^abs(j);
    zeta2 = max(abs(a2(k-j+2*N+1))./wj);
  end
  zeta(k+N+1) = max(zeta1, zeta2);
end
% for k=-N:-1
%   j = (k-2*N:-N-1)'; wj = nu.^abs(j);
%   zeta(k+N+1) = max(abs(a2(k-j+2*N+1))./wj);
% end
% for k=1:N
%   j = (N+1:k+N)'; wj = nu.^abs(j);
%   zeta(k+N+1) = max(abs(a2(k-j+2*N+1))./wj);
% end
wn = nu^(N+1);
ta_norm = wnorm(ba,nu);
Z1_0 = abs(A(1,1))/wn + abs(A_a0)*zeta;
Z1_1 = wnorm(A_o1,nu)/wn + wnorm(abs(A_a1)*zeta,nu)...
  + 3*abs(epsilon)/abs(mu_k(N+1,bomega))*ta_norm^2;
Z1 = max(Z1_0,Z1_1)

% Z2
k = -N:N;
tA = abs(k).*abs(A_a0);
tB = (k.^2).*abs(A_a0);
tA_norm = wsnorm(tA,nu);
tB_norm = wsnorm(tB,nu);
A_norm = wsnorm(A_a0,nu);
alpha = abs(alpha); beta = abs(beta);
gamma = abs(gamma); epsilon = abs(epsilon);

Z2_20 = tB_norm * ((2+beta*tau^2)*ta_norm + 4*bomega)...
  + tA_norm * (gamma + 2*beta*tau + 1)...
  + A_norm * (6*epsilon*ta_norm);
Z2_30 = 3*tB_norm + 3*epsilon*A_norm;

tA = abs(k).*abs(A_a1);
tB = (k.^2).*abs(A_a1);
tA_bopnorm = bopnorm(tA,(N+1)/abs(mu_k(N+1,bomega)),nu);
tB_bopnorm = bopnorm(tB,1/(bomega^2-(alpha+beta)/((N+1)^2)),nu);
A_bopnorm = bopnorm(A_a1,1/abs(mu_k(N+1,bomega)),nu);

Z2_21 = tB_bopnorm * ((2+beta*tau^2)*ta_norm + 4*bomega)...
  + tA_bopnorm * (gamma + 2*beta*tau + 1)...
  + A_bopnorm * (6*epsilon*ta_norm);
Z2_31 = 3*tB_bopnorm + 3*epsilon*A_bopnorm;


Z2_2 = max(Z2_20, Z2_21)
Z2_3 = max(Z2_30, Z2_31)

Z2 = @(r) Z2_3*r.^2 + Z2_2*r;
%% a root of the radii polynomial
p = @(r) Z2(r).*r - (1-Z1-Z0)*r + Y0;
r = fzero(p,1e-10)
if r<0
  disp('Newton-Kantorovich theorem is not holded.')
end
% figure
fplot(p,[0,0.004],'linewidth',1.6), hold on,
plot([0,0.004], [0,0], 'k-.')
title('The radii polynomial')
xlabel('$r$','interpreter','latex'), ylabel('$p(r)$','interpreter', 'latex')
hold off
