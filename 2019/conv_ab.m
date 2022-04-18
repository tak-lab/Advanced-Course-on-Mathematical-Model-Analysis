function ab = conv_ab(a,b)
N = (length(a) - 1)/2;
M = N+1; % M = N+1

p = 2;

N_pad = (p-1)*M;

ta1 = [zeros(N_pad,1);a;zeros(N_pad,1)];% 1. Padding zeros
ta2 = [zeros(N_pad,1);b;zeros(N_pad,1)];% 1. Padding zeros

tb1 = ifft(ifftshift(ta1));% 2. IFFT of ta
tb2 = ifft(ifftshift(ta2));% 2. IFFT of ta

tab = tb1.*tb2;% 3. tb*^tb

ab = fftshift(fft(tab))*(2*p*M-1)^(p-1);% 4. FFT of tb2
ab = ab(N_pad+1:end-N_pad);
