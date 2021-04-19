function plot_solution(u,index)
% index = 1: profile of solution
%         2: Fourier mode
%         3: phase profile
L = 2*pi/real(u(1));
a = u(2:end);
m = length(u)/2;
m_pad = 1000;
a_pad = [zeros(m_pad,1);a;zeros(m_pad,1)];
N = m-1;
N_pad = m+m_pad-1;
k = (-N_pad:N_pad)';
dx = L/(2*N_pad-1);
x = dx*(0:2*N_pad);
if index==1
% Plot profile:
% figure
    plot(x,real((2*N_pad+1)*ifft(ifftshift(a_pad))),'linewidth',1.6)
% hold on
xlabel('$t$','interpreter','latex'), ylabel('$x(t)$','interpreter', 'latex')
% axis([0,2,0,5])
% title('Real part')
elseif index==2
   semilogy(-N:N,abs(a),'linewidth',1.6)
elseif index==3
  % Plot phase:
  plot(real((2*N_pad+1)*ifft(ifftshift(a_pad))),real((2*N_pad+1)*ifft(ifftshift(a_pad.*(1i*k)))),'linewidth',1.6)
% hold on
xlabel('$x$','interpreter','latex'), ylabel('$\dot{x}$','interpreter', 'latex')
% axis([0,2,0,5])
% title('Real part')
end