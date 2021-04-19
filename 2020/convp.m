function cc2_cc = convp(cc,p)
M = length(cc);
ac = [flipud(0.5*cc(2:end));cc(1);0.5*cc(2:end)];% Fourier coefficient
N = (p-1)*M;

ta = [zeros(N,1);ac;zeros(N,1)];% 1. Padding zeros
tb = ifft(ifftshift(ta));% 2. IFFT of ta
tb2 = tb.^p;% 3. tb*^tb
cc2 = real(fftshift(fft(tb2)))*(2*p*M-1)^(p-1);% 4. FFT of tb2
cc2_cc = [cc2(N+M);2*cc2(N+M+1:N+M+p*(M-1))];% Take Chevshev coefficient
