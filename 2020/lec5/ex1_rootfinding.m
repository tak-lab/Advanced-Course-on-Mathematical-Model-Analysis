%% Problem f(x) = 0
f = @(x) x^2 - 2;
DF = @(x) 2*x;

tol = 1e-12;

%% Newton iteration
x0 = 1.0;
x = x0;
while 1
    dx = -f(x)/DF(x);
    x = x + dx;
    if abs(dx) < tol
        break
    end
end
disp(x)

%% A posteriori validation 
Y0 = 1/2/abs(x)*abs(f(x))
Z0 = abs(1-(1/x)*x)
Z1 = 0
Z2 = 1/abs(x)

% p(r) = Z2*r^2 - (1-Z1-Z0)*r + Y0
a = Z2;
b = -(1-Z1-Z0);
c = Y0;
if b^2-4*a*c < 0
    disp('error: cannot find root of radii polynomial')
else
    r = (-b - sqrt(b^2-4*a*c))/2/a
%     r = 2*c/(-b + sqrt(b^2-4*a*c));
%     r = (-b + sqrt(b^2-4*a*c))/2/a
%     r = 2*c/(-b - sqrt(b^2-4*a*c))
end
% disp(r)