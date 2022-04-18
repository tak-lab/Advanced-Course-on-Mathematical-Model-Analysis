N=100;
n=N-1;
x = (0:N-1)'/n*2*pi;
u0=exp(-100*(x-1).^2);
plot(x,u0)
[t,y]=ode45(@advec,[0,1],u0);
mesh(x,t,y)

xlabel 'x'
ylabel 't'
zlabel 'u'