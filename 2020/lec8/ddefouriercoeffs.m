function FourierCoeffs = ddefouriercoeffs(f,I,N,n)
a = I(1); b = I(2);
% x_j: equidistance node points
h = (b-a)/(2*N+1);
j = 0:2*N;
x_j = a + j*h;
% f_j: function values on node points
f_j = deval(f,x_j,n);
% plot(f.x, f.y, LW, lw), hold on
% plot(x_j,f_j,'ro')
FourierCoeffs = (fftshift(fft(f_j)))/(2*N+1);