
N=8; %period of DFT
s = randn(1,8);
h = randn(1,2);

cir_s_h = cconv(h,s,N) % circular convolution of h and s with period N

Ncp = 2;
s_cp = [s(end-Ncp+1:end) s]; %copy last Ncp symbols from s and prefix it.
lin_s_h = conv(h,s_cp)

H = fft(h,N); %frequency response of CIR
S = fft(s,N); %frequency response of OFDM signal (non CP)
 
r1 = ifft(S.*H) %IFFT of product of individual DFTs