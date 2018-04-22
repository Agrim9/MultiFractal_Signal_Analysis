% This is a function implementation of wavelets_opnheim.m
function [gam_est, sig_x, sig_w, sig_reg, var_w] = oppnheim_func(gamma, N, snr, wavelet)
M = log2(N);
cn = dsp.ColoredNoise(gamma, N,1); 
rng default; %for repeatability
x = step(cn); % generating 1/f process

r = awgn(x, snr,'measured'); %adding white gaussian noise

[r_mn, l1] = wavedec(r, M, wavelet);
r_m = detcoef(r_mn, l1, 1:M);

power_x = rms(x)^2;
var_w = power_x/db2pow(snr);

[x_mn, l] = wavedec(x,M, wavelet); %finding wavelet transform of x
x_m = detcoef(x_mn, l, 1:M);

varx = zeros(M,1);
for i = 1:M
    varx(M+1-i)=var(x_m{1,i}); % variance vector inverted for correct seq. of m
end

% estimating sigma and gamma using linear regression

B = [ones(ceil(0.5*M),1) ((floor(0.5*M)+1):M)']\log2(varx((floor(0.5*M)+1):M)); 
sig_reg = 2^B(1);

% finding variance at diff scales of noisy signal
varr = zeros(M,1);
for i = 1:M
    varr(M+1-i)=var(r_m{1,i});
end

% starting E-M algorithm

[beta, sig_x, sig_w] = EM_estimate(1, 1, 1, l1, varr); 
gam_est = log2(beta);

end