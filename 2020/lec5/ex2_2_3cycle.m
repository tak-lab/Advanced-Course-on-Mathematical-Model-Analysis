clear
%% Newton iteration
tol = 1e-12;
p = 1;
% x0 = [3.83;1;-1;1];
x0 = [3;1;-1;1];
% eta = 1.6429;
eta = 1.6311;
% eta = 1.75;
% eta = 1.65;
x = x0;
% disp(['Before iteration: ', num2str(norm(F(x,eta)))])
step = 1;
while 1
    dx = -DF(x)\F(x,eta);
    x = x + dx;
%     disp(['After ',num2str(step),' iteration: ', num2str(norm(F(x,eta)))])
    step = step + 1;
    if or(norm(dx,p) < tol, step>500)
        break
    end
end
disp(x)

%% A posteriori validation
A = inv(DF(x));
Y0 = norm(A*F(x,eta),p)
Z0 = norm(eye(size(A))-A*DF(x),p)
Z1 = 0
lam = x(1); x(1) = [];
B = @(r) [0,0,0,0;...
  (1+(2*x(3)+r)), 0,0,(1+(x(3)+r))+(lam+r);...
  (1+(2*x(1)+r)),(1+(x(1)+r))+(lam+r),0,0;...
  (1+(2*x(2)+r)),0, (1+(x(2)+r))+(lam+r),0];
Z2 = @(r) norm(A*B(r),p);
% 
p = @(r) Z2(r)*r^2 - (1-Z1-Z0)*r + Y0;
r = fzero(p,0.01)
% disp(r)

fplot(p,[0,0.06],'linewidth',1.6)
title('The radii polynomial')
xlabel('$r$','interpreter','latex'), ylabel('$p(r)$','interpreter', 'latex')
hold on
plot([0,0.06],[0,0],'r--')

function F=F(x,eta)
lam = x(1); x(1) = [];
F = [sum(x)-eta;...
  x(1) - lam*x(3)*(1-x(3));...
  x(2) - lam*x(1)*(1-x(1));...
  x(3) - lam*x(2)*(1-x(2))];
end

function DF = DF(x)
lam = x(1); x(1) = [];
DF = [0, 1, 1, 1;...
  -x(3)*(1-x(3)), 1, 0, -lam*(1-2*x(3));...
  -x(1)*(1-x(1)), -lam*(1-2*x(1)), 1, 0;...
  -x(2)*(1-x(2)), 0, -lam*(1-2*x(2)), 1];
end
