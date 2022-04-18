tend=10;
for i=linspace(.05,3,30)
  [t1,y1] = ode45(@(t,y)y*(1-y),[0,tend],i);
  plot(t1,y1)
  hold on
end