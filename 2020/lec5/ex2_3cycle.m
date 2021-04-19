clear
%% Problem f(x) = 0
F = @(x,lam) [x(1) - lam*x(3)*(1-x(3));x(2) - lam*x(1)*(1-x(1));x(3) - lam*x(2)*(1-x(2))];
DF = @(x,lam) [1, 0, -lam*(1-2*x(3)); -lam*(1-2*x(1)), 1, 0; 0, -lam*(1-2*x(2)), 1];

tol = 1e-12;
lam = 3.82843; % stable 3 cycle in the 3-cycle window
% lam = 4; % unstable 3 cycle 
p = 1;

%% Newton iteration
x0 = [1;-1;1];
x = x0;
disp(['Before iteration: ', num2str(norm(F(x,lam)))])
step = 1;
while 1
    dx = -DF(x,lam)\F(x,lam);
    x = x + dx;
    disp(['After ',num2str(step),' iteration: ', num2str(norm(F(x,lam)))])
    step = step + 1;
    if norm(dx,p) < tol
        break
    end
end
disp(x)

%% A posteriori validation
A = inv(DF(x,lam));
Y0 = norm(A*F(x,lam),p)
Z0 = norm(eye(size(A))-A*DF(x,lam),p)
Z1 = 0
Z2 = 2*abs(lam)*norm(A,p)

% p(r) = Z2*r^2 - (1-Z1-Z0)*r + Y0
a = Z2;
b = -(1-Z1-Z0);
c = Y0;
if b^2-4*a*c < 0
    disp('error: cannot find root of radii polynomial')
else
    r = (-b - sqrt(b^2-4*a*c))/2/a;
%     r = 2*a*c/(-b + sqrt(b^2-4*a*c));
end
disp(r)
