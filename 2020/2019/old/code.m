% MATLABの接続方法
% https://www.u.tsukuba.ac.jp/remote/
% https://www.u.tsukuba.ac.jp/remote/#view2win

% コード一覧
% 2018/04/10
a=-1;f=@(t,x)a*x;
[t,y]=ode45(f,[0,3],1);
LW = 'linewidth'; lw = 1.6;
plot(t,y,LW,lw)

% 2018/04/17
beta=2;nu=1;f=@(t,y,beta,nu) [-beta*y(1)*y(2);(beta*y(1)-nu)*y(2)];
[t,y]=ode45(@(t,y)f(t,y,beta,nu),[0,30],[1-0.0001,0.0001]);
plot(t,[y,1-y(:,1)-y(:,2)],LW,lw)
figure
plot(y(:,1),y(:,2),LW,lw)
for beta=[1.3, 1.6, 2, 2.8]
f=@(t,y,beta,nu) [-beta*y(1)*y(2);(beta*y(1)-nu)*y(2)];
[t,y]=ode45(@(t,y)f(t,y,beta,nu),[0,30],[1-0.0001,0.0001]);
plot(t,y(:,2),LW,lw), hold on
end
for i=1:-0.1:0.1
for init=[0.1:0.1:0.9]
[t,y]=ode45(@(t,y)f(t,y,beta,nu,mu),[0,30],[i-init,init]);
plot(y(:,1),y(:,2)), hold on
end
end
axis([0,1,0,1])
plot(1/2,1/4,'o',LW,lw)

% 2018/04/24
f=@(t,y,a,b,c,d) [a*y(1)-b*y(1)*y(2);-c*y(2)+d*y(1)*y(2)];
[t,y]=ode45(@(t,y) f(t,y,1,0.01,1,0.02),[0,30],[1,1]);
LW = 'linewidth'; lw = 1.6;
plot(t,y,LW,lw)
[t,y]=ode45(@(t,y) f(t,y,1,1,0.2,1),[0,60],[1,1]);
plot(t,y,LW,lw)
plot(y(:,1),y(:,2),LW,lw), hold on
plot(0.2,1,'ro')
plot([0,1.2],[1,1],'--')
plot([0.2,0.2],[0,2.5],'--')
% Chebfun codes
tol=1e-20;opts = odeset('abstol',tol,'reltol',tol);
y=chebfun.ode45(@(t,y) f(t,y,a,b,c,d),[0,60],[1,1],opts);
plot(y(:,1),y(:,2),LW,lw)
arrowplot(y(:,1),y(:,2),LW,lw)
plot(y(:,1),LW,lw)
[ym,m]=max(y(:,1),'local');
hold on
plot(m,ym,'ro')
diff(m)
N=chebop(0,60);
N.lbc=[1;1];
N.op = @(t,u,v)[diff(u)-u+u*v;diff(v)+0.2*v-u*v];
quiver(N,[0 1.5 0 3])
hold on
plot([0,1.2],[1,1],'--')
plot([0.2,0.2],[0,2.5],'--')
f=@(t,y,a,b,c,d) [a*y(1)-b*y(1)*y(2);-c*y(2)+d*y(1)*y(2)];
a=1;b=1;c=0.2;d=1;
for s = .4:.2:1.4
y=chebfun.ode45(@(t,y) f(t,y,a,b,c,d),[0,60],[s,1],opts);
arrowplot(y(:,1),y(:,2),LW,lw), hold on,
end
plot([0.2,0.2],[0,],'--')
plot([0.2,0.2],[0,3],'--')
plot([0,1.4],[1,1],'--')

% 2018/05/01
f=@(t,y,sig,r,b) [sig*(y(2)-y(1));-y(2)+y(1)*(r-y(3));y(1)*y(2)-b*y(3)];
sig=1;r=0.5;b=1;
[t,y]=ode45(@(t,y) f(t,y,sig,r,b),[0,30],[1,1,1]);
LW = 'linewidth'; lw = 1.6;
plot(t,y,LW,lw)
plot3(y(:,1),y(:,2),y(:,3),LW,lw)
for i=100:100:1000
[t,y]=ode45(@(t,y) f(t,y,sig,r,b),[0,30],[i*randn,i*randn,i*randn]);plot3(y(:,1),y(:,2),y(:,3),LW,lw)
hold on
pause
end
sig=10;r=20;b=8/3;
rs=sig*(sig+b+3)/(sig-b-1);
[t,y]=ode45(@(t,y) f(t,y,sig,r,b),[0,30],[1,1,1]);
plot(t,y,LW,lw)
plot3(y(:,1),y(:,2),y(:,3),LW,lw)
plot3(-sqrt(b*(r-1)),-sqrt(b*(r-1)),r-1,'ro',LW,lw)
sig=10;r=28;b=8/3;
[t,y]=ode45(@(t,y) f(t,y,sig,r,b),[0,30],[1,1,1]);
plot(t,y,LW,lw)
plot3(y(:,1),y(:,2),y(:,3),LW,lw)
hold on
plot3(1,1,1,'ro',LW,lw)
plot3(y(end,1),y(end,2),y(end,3),'vr',LW,lw)
[t,y]=ode45(@(t,y) f(t,y,sig,r,b),[0,30],[-1,-1,1]);
plot3(y(:,1),y(:,2),y(:,3),LW,lw)
[t,y]=ode45(@(t,y) f(t,y,sig,r,b),[0,30],[-1,-1,-1]);
plot3(y(:,1),y(:,2),y(:,3),LW,lw)
[t,y]=ode45(@(t,y) f(t,y,sig,r,b),[0,30],[-15,-15,20]);
plot3(y(:,1),y(:,2),y(:,3),LW,lw), hold on
[t,y]=ode45(@(t,y) f(t,y,sig,r,b),[0,30],[-15,-15,20.00001]);
plot3(y(:,1),y(:,2),y(:,3),LW,lw)
[t,y]=ode45(@(t,y) f(t,y,sig,r,b),[0,30],[-15,-15,20]);
[t1,y1]=ode45(@(t,y) f(t,y,sig,r,b),[0,30],[-15,-15,20.00001]);
plot(t,y(:,2),LW,lw),hold on,plot(t1,y1(:,2),LW,lw)
tol=1e-20;opts = odeset('abstol',tol,'reltol',tol);
[t2,y2]=ode45(@(t,y) f(t,y,sig,r,b),[0,30],[-15,-15,20],opts);
plot(t,y(:,2),LW,lw),hold on,plot(t2,y2(:,2),LW,lw)
N=chebop(0,30);
N.op=@(t,u,v,w)[diff(u)-10*(v-u);diff(v)-u*(28-w)+v;diff(w)-u*v+(8/3)*w];
N.lbc=@(u,v,w)[u+15; v+15; w-20];
[u,v,w] = N\0;
plot(v,LW,lw)
fun = @(t,u) [10*(u(2)-u(1)); 28*u(1)-u(2)-u(1)*u(3); u(1)*u(2)-(8/3)*u(3)];
u = chebfun.ode113(fun,[0,30],[-15 -15 20], opts);
plot(u(:,2),LW,lw)
% 精度は\よりode45やchebfun.ode113が勝利する模様

% 2018/05/08
n=20;x=(0:(1/n):1)';x(1)=[];x(end)=[];
f=@(x) sin(pi*x)+sin(3*pi*x)+sin(12*pi*x);
[t,y]=ode45(@(t,y) func(t,y),[0,1],f(x));
mesh(x,t',y)
for k=1:100
  plot([0;x;1],[0,y(k,:),0],LW,lw)
  axis([0,1,-1,3])
  pause(0.1)
end
[t,y]=ode45(@(t,y) func(t,y),[0,1],100*rand(size(x)));
[t,y]=ode45(@(t,y) func(t,y),[0,1],f(x));
u=chebfun.ode45(@(t,y) func(t,y),[0,1],f(x));
tspan=chebpts(length(u),[0,1]);
mesh(x,tspan,u(tspan))
hold on
mesh(x,t',y)

% 2018/05/15
func1 = @(t,y) fujita_neumann(t,y,3);
f=@(x) sin(pi*x);
n=20;x=(0:(1/n):1)';x(1)=[];x(end)=[];
[t,y]=ode45(@(t,y) func1(t,y),[0,1],f(x));
mesh(x,t',y)
LW = 'linewidth'; lw = 1.6;
plot(x,y(1,:),LW,lw)
[t,y]=ode45(@(t,y) func1(t,y),[0,1.2],f(x));
mesh(x,t',y)
set(gca,'ZScale','log')
u=chebfun.ode45(@(t,y) func1(t,y),[0,1.2],f(x));
u=chebfun.ode45(@(t,y) func1(t,y),[0,1.11214],f(x));
tspan=chebpts(length(u),[0,1.11214]);
mesh(x,tspan,u(tspan))
set(gca,'ZScale','log')
figure
mesh(x,t',y)
set(gca,'ZScale','log')
f=@(x) sin(pi*x)+sin(3*pi*x)+sin(12*pi*x);
plot(x,f(x),LW,lw)
[t,y]=ode45(@(t,y) func1(t,y),[0,1.2],f(x));
mesh(x,t',y)
set(gca,'ZScale','log')
mesh(x,t',y)
plot(x,y(1,:),LW,lw)
mesh(x,t',y)
n=100;x=(0:(1/n):1)';x(1)=[];x(end)=[];
plot(x,f(x),LW,lw)
[t,y]=ode45(@(t,y) func1(t,y),[0,1.2],f(x));
mesh(x,t',y)
set(gca,'ZScale','log')
[t,y]=ode45(@(t,y) func1(t,y),[0,5],rand(size(x)));
mesh(x,t',y)
set(gca,'ZScale','log')
plot(x,y(1,:),LW,lw)

% 2018/6/5
xy = @(t) 1+t-4*t.^2+5*t.^3-3*t.^4;
LW = 'linewidth'; lw = 1.6;
fplot(xy,[0,0.1],LW.lw)
hold on
x_1 = @(t) 1+t-4*t.^2;x_2 = @(t) 1+t-3.5*t.^2;
fplot(x_1,[0,0.1],LW,lw)
fplot(x_2,[0,0.1],LW,lw)
clf
fplot(xy,[0,0.5],LW,lw),hold on
fplot(x_1,[0,0.5],LW,lw)
fplot(x_2,[0,0.5],LW,lw)