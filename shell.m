%shell to obtain plots for varying gamma and N and snr

gamma = [0.75 1 1.25];
snr = 0:5:50;
N = [2^8 2^10 2^12 2^14];
wavelet = ['db2db3db4db5'];
K = size(snr,2);
L = size(wavelet,2)/3;

% gam_est = zeros(K,1);
% sig_x = zeros(K,1);
% sig_w = zeros(K,1);
% sig_act = zeros(K,1);
% var_w = zeros(K,1);

gam_est = zeros(K,L);
sig_x = zeros(K,L);
sig_w = zeros(K,L);
sig_act = zeros(K,L);
var_w = zeros(K,L);

for j = 1:L
    for i=1:K
        [gam_est(i,j), sig_x(i,j), sig_w(i,j), sig_act(i,j), var_w(i,j)] = oppnheim_func(1, 8192, snr(i), wavelet((3*j-2):3*j));
    end
end
%     for i=1:K
%         [gam_est(i), sig_x(i), sig_w(i), sig_act(i), var_w(i)] = oppnheim_func(gamma(2), N(3), snr(i), 'db5');
%     end

% gamma_err = abs(gam_est - ones(K,1)*gamma);

gamma_err = abs(gam_est - ones(K,L));
% 
% plot((snr)',gamma_err(:,1), 'marker', 'o');
% grid on
% hold on
% plot((snr)',gamma_err(:,2), 'marker', '^');
% plot((snr)',gamma_err(:,3),'marker', 'd');
% hold off
% ylim([0 0.4])
% title('Error in \gamma with SNR')
% legend('\gamma = 0.75','\gamma = 1', '\gamma = 1.25')
% xlabel('SNR')
% ylabel('RMS error in \gamma')

plot((snr)',gamma_err(:,1), 'marker', 'o');
grid on
hold on
plot((snr)',gamma_err(:,2), 'marker', '^');
plot((snr)',gamma_err(:,3),'marker', 'd');
plot((snr)',gamma_err(:,4),'marker', 'x');
hold off
ylim([0 0.4])
title('Error in \gamma with SNR for different wavelets')
legend('db2','db3','db4', 'db5')
xlabel('SNR')
ylabel('RMS error in \gamma')


% semilogx((N)',gamma_err(:,1), 'marker', 'o');
% grid on
% hold on
% semilogx((N)',gamma_err(:,2), 'marker', '^');
% semilogx((N)',gamma_err(:,3),'marker', 'd');
% hold off
% ylim([0 0.4])
% title('Error in \gamma with N')
% legend('\gamma = 0.75','\gamma = 1', '\gamma = 1.25')
% xlabel('Data Length/ Resolution log scale (log(N))')
% ylabel('RMS error in \gamma')