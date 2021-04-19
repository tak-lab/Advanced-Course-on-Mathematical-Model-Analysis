function [ap, c] = convp(a,p)
M = (length(a) +1)/2;
n = (p-1)*M;
ta = [zeros(n,1); a; zeros(n,1)];
tb = ifft(ifftshift(ta));
tc  = fftshift(fft(tb.^p))*(2*p*M-1)^(p-1);
ap = tc(n+1:end-n);
 c = tc(p:end-p+1);