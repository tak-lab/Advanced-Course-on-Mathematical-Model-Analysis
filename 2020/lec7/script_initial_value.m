% Delayed Duffing equation
% % alpha = -0.039;
% alpha = -0.1;
% % beta = -0.4;
% beta = -5;
% gamma = -2;
% epsilon = -10;
% tau = 2.5;

%%
dde = @(t,y,z) [y(2); gamma*y(2) + alpha*y(1) + beta*z(1) + epsilon*y(1).^3];
ddehist = @(t) [0.4;0];
% ddehist = dde23(@(t,y,z) zeros(size(y)),[], [0.5;0], [-tau,0]);
% opts=ddeset('RelTol',1.0e-7,'AbsTol',1.0e-8);
% opts = ddeset('AbsTol',1e-18,'RelTol',1e-18);
% opts = ddeset('RelTol',1e-14);
tspan = [0,400];
delay = tau;
f = dde23(dde,delay,ddehist,tspan);
plot(f.y(1,floor(length(f.x)/2):end).',f.y(2,floor(length(f.x)/2):end).')
% % plot(f.y(1,:),f.y(2,:))

%% approximate a period of solution
a = 200;
app_period = 6;
h=0.1;
f_tmp = deval(f,a+app_period/2:h:a+3*app_period/2);
find_period = abs(f_tmp - deval(f,a));
eta0 = deval(f,a,1);
[~,ind] = min(find_period(1,:));
b = a+app_period/2 + h*(ind-1);

%% transform the data into Fourier coefficient
x = ddefouriercoeffs(f,[a,b],N,1);
x0 = [2*pi/(b-a);x.'];
% figure
plot_solution(x0,1)