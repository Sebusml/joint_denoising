clear all
close all
clc

% load 'data/sensor_data.mat'
% 
% X = parse_sensor_data(sensor_data);
% plot(X(:,1:end))


% load data/temperatures.mat
% 
% i = 1:size(X,2);
% i([2 19]) = [];
% X = X(60:end,:);
% X = X(:,i);
% days = 2;
% n = 1:24*60*2*days;
% X = X(n,:);
% save 'data/temperature_2day.mat' X

load data/temperature_2day.mat
X = X(1:2^floor(log2(size(X,1))),:); % Truncate X to have power of two length
J = size(X,2);
% theta = fft(X);
% Remove the mean value of the signals
% theta(1,:) = zeros(1,J);
% X = real(ifft(theta));

sigma = 0.4;
N = size(X,1);
noise = sigma * randn(N,J);
X_n = X + noise;
T = sigma * sqrt(2*log(N));
T = T * sqrt(N); % The fft is not orthonormal, only orhtogonal. The noise
                 % standard deviation is multiplyed by sqrt(N).
% plot(X_n(:,1)), hold on, plot(X(:,1),'r')
theta = fft(X_n);

thresh_coeff = (abs(theta) > T);
save_coeff = sum(abs(theta)>T,2) > 0;
save_coeff = repmat(save_coeff,1,J);
keep_coeff = (thresh_coeff + save_coeff) > 0;

theta_hat_indep = theta .* thresh_coeff; 
theta_hat_joint = theta .* keep_coeff;

X_hat_indep = real(ifft(theta_hat_indep));
X_hat_joint = real(ifft(theta_hat_joint));

% figure, hold on
% plot(X(:,1),'g')
% plot(X_hat_indep(:,1),'r')
% plot(X_hat_joint(:,1),'b')
% legend('Original', 'Indep', 'Joint')

risk_indep = norm(X(:,1) - X_hat_indep(:,1))^2
risk_joint = norm(X(:,1) - X_hat_joint(:,1))^2


%  Plots
ns  = 1:N;
lw = 1; % Linewidth
yl = [16 26];
% m = max(abs(X_n(:,1))) + 0.4;
figure(1)
plot(ns,X(:,1),'b','LineWidth',lw)
xlabel('$k$')
ylabel('Temperature $[^\circ C]$')
title('Original signal ($j=1$)')
xlim([1 N])
ylim(yl)

figure(2)
plot(X_n(:,1),'linewidth',lw)
xlabel('$k$')
ylabel('Temperature $[^\circ C]$')
title(['Noisy signal $\sigma_w$ = ' num2str(sigma)])
xlim([1 N])
ylim(yl)

figure(3)
xj = X_hat_joint(:,1);
xi = X_hat_indep(:,1);
plot(ns,X(:,1),'k:','linewidth',lw)
hold on
plot(ns, xi,'b','linewidth',lw)
plot(ns, xj,'r','linewidth',lw)
legend('Original','Independent','Joint')
ylabel('Temperature $[^\circ C]$')
xlabel('$k$')
title('Denoised signals')
xlim([1 N])
ylim(yl)

%% Save data
save_data = true;
if save_data,
    save ('jd_temperatures', '-V7')
end
% plotpdftex(1,'exp_1',[1 1])
% plotpdftex(2,'exp_2',[1 1])
% plotpdftex(3,'exp_3',[1 1])
