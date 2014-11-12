clear all
close all
clc

% Analyze the improvement equation

%%
N = 1024;
sparsity_ratio = 0.05;
S = round(sparsity_ratio * N);
sigma_x = 1;
sigma_w = 0.4;%0.4;
T = sigma_w * sqrt(2*log(N));
J = 2:1:100;

imp = improvement(T, J, sigma_x, sigma_w, S, N);
figure(1)
plot(J, imp, 'LineWidth',2)
xlabel('Number of signals')
ylabel('Expected improvement')
s = sprintf('$N=%d$, $S=%d$, $\\sigma_x=%.1f$, $\\sigma_w=%.1f$, $T=%.1f$', ...
        N, S, sigma_x, sigma_w, T);
title(s)

%plotpdftex(1,'improvement',[1.4 1])

%%
N = 1024;
sparsity_ratio = 0.05;
S = round(sparsity_ratio * N);
sigma_x = 1;
sigma_w = linspace(0.1,1.2,50);
T = sigma_w * sqrt(2*log(N));
J = 100;

imp = improvement(T, J, sigma_x, sigma_w, S, N);
figure(2)
plot(sigma_w, imp, 'LineWidth',2)
xlabel('Noise variance')
ylabel('Expected improvement')
grid
% s = sprintf('$N=%d$, $S=%d$, $\\sigma_x=%.1f$, $\\sigma_w=%.1f$, $T=%.1f$', ...
%         N, S, sigma_x, sigma_w, T);
% title(s)

N = 1024;
sparsity_ratio = 0.05;
S = round(sparsity_ratio * N);
sigma_x = linspace(1,10,50);
sigma_w = 0.4;
T = sigma_w * sqrt(2*log(N));
J = 10;

imp = improvement(T, J, sigma_x, sigma_w, S, N);
figure(3)
plot(sigma_x, imp, 'LineWidth',2)
xlabel('Signal variance')
ylabel('Expected improvement')
grid


% T = T * sqrt(log10(J));
% for i = 1:length(J),
%     imp(i) = improvement(T(i), J(i), sigma_x, sigma_w, S, N);
% end
% % subplot(212)
% hold on
% plot(J, imp,'r')

% figure
% plot(T)

