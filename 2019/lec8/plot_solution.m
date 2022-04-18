function plot_solution(u,index)
% index = 1: profile of solution
%         2: Fourier mode
%         3: profile of solution on phase space

m = length(u); % m = 4*N+3
N = (m-3)/4;

L = 2*pi/real(u(1));
a = u(2:2*N+2);
b = u(2*N+3:end);

dx = L/(2*N-1);
x = dx*(0:2*N);
if index==1
% Plot profile:
% figure
    subplot(1,2,1)
    plot(x,real((2*N+1)*ifft(ifftshift(a))),'linewidth',1.6)
xlabel('$t$','interpreter','latex'), ylabel('$x$','interpreter', 'latex')
    subplot(1,2,2)
    plot(x,real((2*N+1)*ifft(ifftshift(b))),'linewidth',1.6)
xlabel('$t$','interpreter','latex'), ylabel('$y$','interpreter', 'latex')
% hold on
% xlabel('$x$','interpreter','latex'), ylabel('$t$','interpreter', 'latex'), zlabel('$\mathrm{Re}(u)$','interpreter', 'latex')
% title('Real part')
elseif index==2
   subplot(1,2,1)
   semilogy(-N:N,abs(a),'linewidth',1.6)
xlabel('$k$','interpreter','latex'), ylabel('$|a_k|$','interpreter', 'latex')
   subplot(1,2,2)
   semilogy(-N:N,abs(b),'linewidth',1.6)
xlabel('$k$','interpreter','latex'), ylabel('$|b_k|$','interpreter', 'latex')
elseif index==3
    plot(real((2*N+1)*ifft(ifftshift(a))),real((2*N+1)*ifft(ifftshift(b))),'linewidth',1.6)
xlabel('$x$','interpreter','latex'), ylabel('$y$','interpreter', 'latex')
end