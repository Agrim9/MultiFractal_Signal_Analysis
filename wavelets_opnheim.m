% This code is an implementation of the paper by Wornell and Oppenheim
gamma = 1;
N = 65536;
M = log2(N);
snr= 0;

cn = dsp.ColoredNoise(gamma, N,1); 
rng default; %for repeatability
x = step(cn); % generating 1/f process
% calculating parameters of signal
variance_x = var(x);
power_x = rms(x)^2;
variance_w = power_x/db2pow(snr);

r = awgn(x, snr,'measured'); %adding white gaussian noise
var_r = var(r);

[x_mn, l] = wavedec(x,M,'db5'); %finding wavelet transform of x
x_m = detcoef(x_mn, l, 1:M);

varx = zeros(M,1);
for i = 1:M
    varx(M+1-i)=var(x_m{1,i}); % variance vector inverted for correct seq. of m
end

% estimating sigma and gamma using linear regression

B = [ones(ceil(0.5*M),1) ((floor(0.5*M)+1):M)']\log2(varx((floor(0.5*M)+1):M)); 
gam_reg = -B(2);
sig_reg = 2^B(1);

% finding wavelet transform of noisy signal
[r_mn, l1] = wavedec(r, M,'db5');
r_m = detcoef(r_mn, l1, 1:M);

% finding variance at diff scales of noisy signal
varr = zeros(M,1);
for i = 1:M
    varr(M+1-i)=var(r_m{1,i});
end
% starting E-M algorithm

[beta, sig_x, sig_w] = EM_estimate(1, 1, 1, l, varr); 
gam_est = log2(beta);

% estimates from direct calculation step of EM algorithm (sig_w=0)
Cm = (1:M)'./sum((1:M)'.*l(2:M+1)) - ones(M,1)./sum(l(2:M+1));
coeff = Cm.*l(2:M+1).*varx ;
rt = roots(flipud(coeff));
beta_est_1 = rt(imag(rt)==0 & real(rt)>0);
sig_est_1 = sum(l(2:M+1).*varx.*beta_est_1.^((1:M)'))/sum(l(2:M+1));
gam_est_1 = log2(beta_est_1);

%reconstruction
x_1 = waverec(x_mn, l,'db5');
    
x_mn_hat = zeros(size(x_mn));
x_mn_hat(1:l(1))= x_mn(1:l(1));

for i=1:M 
    x_mn_hat(l(i)+1:l(i+1)) = (sig_x*beta^-i/(sig_w + sig_x*beta^(-i)))*r_mn(l(i)+1:l(i+1));
end

x_hat = waverec(x_mn_hat, l, 'db5');

%%%%%%%%----PLOTS------%%%%%%%%%%%%

% Regression plots
figure(1)
plot(log2(varx))
hold on 
grid on
plot(B(1)+B(2)*(1:M))
hold off
title('Estimaton of signal variance')
xlabel('scale (m)')
ylabel('log_2(x^m_n)')


%reconstruction plots
figure(2)
title('Signal Reconstruction')
subplot(3,1,1)
plot(x(1:1000))
xlim([0 1000])
ylim([-10 10])
grid on;
xlabel('t')
ylabel('fractal signal x(t)')
subplot(3,1,2)
plot(r(1:1000))
xlim([0 1000])
ylim([-10 10])
grid on;
xlabel('t')
ylabel('signal in noise r(t)')
subplot(3,1,3)
plot(x_hat(1:1000))
xlim([0 1000])
ylim([-10 10])
grid on;
xlabel('t')
ylabel('estimated signal x^{^}(t)')

% subplot(4,1,4)
% plot(x_1(1:1000))
% xlim([0 1000])
% % ylim([-20 20])
% grid on;

% PSD Plots
figure(3)
Fs = 1;
[Pxx,F] = pwelch(x,hamming(128),[],[],Fs,'psd');
PSDPink = 1./F(2:end);
plot(log2(F(2:end)),10*log10(Pxx(2:end)))
hold on
plot(log2(F(2:end)),10*log10(PSDPink),'r','linewidth',2)
xlabel('log_2(Hz)')
ylabel('dB')
title('1/f signal Power Spectrum for \gamma = 1')
grid on
legend('PSD estimate','Theoretical ')
hold off