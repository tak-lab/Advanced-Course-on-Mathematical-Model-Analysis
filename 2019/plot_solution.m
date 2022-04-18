function plot_solution(u,index)
% index = 1: profile of solution
%         2: Fourier mode
L = 2*pi/real(u(1));
a = u(2:end);
m = length(u)/2;
N = m-1;
dx = L/(2*N-1);
x = dx*(0:2*N);
if index==1
% Plot profile:
% figure
    plot(x,real((2*N+1)*ifft(ifftshift(a))),'linewidth',1.6)
% hold on
xlabel('$t$','interpreter','latex'), ylabel('$y(t)$','interpreter', 'latex')
% title('Real part')
elseif index==2
   semilogy(-N:N,abs(a),'linewidth',1.6)
end