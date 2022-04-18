n = 32;
k = 0:n;
x = k/n;
init = 50*(1.-cos(2*pi*x));
[t,y] = ode45('heat',[0,3],init(2:n));

m = size(t,1);
for i=1:m-1
  u = [0,y(i,:),0];
  plot(x,u)
  ylim([0 100])
  pause(t(i+1)-t(i))
end