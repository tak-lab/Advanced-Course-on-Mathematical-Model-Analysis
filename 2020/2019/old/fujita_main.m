n = 32;
k = 0:n;
x = k/n;
init = 10*(1.-cos(5*pi*x));
[t,y] = ode45('fujita_neu',[0,3],init(2:n));

m = size(t,1);
for i=1:m-1
  u = [y(i,1),y(i,:),y(i,end)];
  plot(x,u)
  ylim([-50 50])
  xlim([0,1])
  pause(t(i+1)-t(i))
end