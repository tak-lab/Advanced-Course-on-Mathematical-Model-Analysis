tol = 5e-10;
% alpha = -0.051;
% beta = -0.4;
alpha = -0.1;
beta = -5;
gamma = -2;
epsilon = -10;
tau = 2.5;
% for tau = 5:0.1:8
% initial value is needed
N = 200;
script_initial_value
eta0 = real(sum(x0(2:end)));
x = x0;
% eta0 = 0.07;
for tau=2.5:0.1:8
% 
% map F
F = F_fourier(x,alpha,beta,gamma,epsilon,tau,eta0);
% 
display(['Before step #1, ||F||_1 = ',num2str(norm(F,1))])
itercount = 0;

while (itercount<=100) && (norm(F,1) > tol)
    DF = DF_fourier(x,alpha,beta,gamma,epsilon,tau);
%     DF=finite_diff_DF(x,alpha,beta,gamma,epsilon,tau,eta0);
    x = x - DF\F;
    F = F_fourier(x,alpha,beta,gamma,epsilon,tau,eta0);
    display(['After step # ',num2str(itercount+1),', ||F||_1 = ',num2str(norm(F,1))])
    itercount = itercount+1;
end

if  2*pi/real(x(1))<0
  disp('fail')
end
plot_solution(x,3),pause(0.1),tau
end
% end