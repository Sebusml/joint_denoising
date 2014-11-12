clear all
close all
clc

N = 2048;
J = 50;
sigma = 0.01;

sparsity_rate = 0.05:0.05:0.8;
original = zeros(length(sparsity_rate),1);
indep = zeros(length(sparsity_rate),1);
joint = zeros(length(sparsity_rate),1);
for i = 1:length(sparsity_rate),
    K = round(sparsity_rate(i)*N);
    [SNR_original, SNR_indep, SNR_joint, risk_indep, risk_joint] = ...
        joint_denoising_time_2(J,N,K,sigma);
    original(i) = 20*log10(mean(SNR_original));
    indep(i) = 20*log10(mean(SNR_indep));
    joint(i) = 20*log10(mean(SNR_joint));
end
subplot(311), hold on
plot(sparsity_rate,original,'g')
plot(sparsity_rate,indep,'b')
plot(sparsity_rate,joint,'r')
legend('Original', 'Independent','Joint')
xlabel('Sparsity Rate')
ylabel('SNR[dB]')
s = sprintf('N = %d, J = %d, sigma = %f', N, J, sigma);
title(s)

J = 50;
K = round(0.1*N);

sigma = 0.005:0.001:0.02;
original = zeros(length(sigma),1);
indep = zeros(length(sigma),1);
joint = zeros(length(sigma),1);
for i = 1:length(sigma),
    [SNR_original, SNR_indep, SNR_joint] = ...
        joint_denoising_time_2(J,N,K,sigma(i));
    original(i) = 20*log10(mean(SNR_original));
    indep(i) = 20*log10(mean(SNR_indep));
    joint(i) = 20*log10(mean(SNR_joint));
end
subplot(312), hold on
plot(sigma,original,'g')
plot(sigma,indep,'b')
plot(sigma,joint,'r')
legend('Original', 'Independent','Joint')
xlabel('Noise sigma')
ylabel('SNR[dB]')
s = sprintf('N = %d, J = %d, K = %f', N, J, K);
title(s)

K = round(0.1*N);
sigma = 0.01;

J = 2:1:20;
original = zeros(length(J),1);
indep = zeros(length(J),1);
joint = zeros(length(J),1);
for i = 1:length(J),
    [SNR_original, SNR_indep, SNR_joint] = ...
        joint_denoising_time_2(J(i),N,K,sigma);
    original(i) = 20*log10(mean(SNR_original));
    indep(i) = 20*log10(mean(SNR_indep));
    joint(i) = 20*log10(mean(SNR_joint));
end
subplot(313), hold on
plot(J,original,'g')
plot(J,indep,'b')
plot(J,joint,'r')
legend('Original', 'Independent','Joint')
xlabel('Number of signals')
ylabel('SNR[dB]')
s = sprintf('N = %d, K = %d, sigma = %f', N, K, sigma);
title(s)