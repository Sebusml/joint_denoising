
clear all
close all
clc

J = 200;
N = 2^10;
K = 2^6;
X = make_jsm2_signals(J,N,'time',K,'normal');
sigma = 0.04;
X_n = X + sigma*randn(N,J);
T = sigma * sqrt(2*log(N));


thresh_coeff = (abs(X_n) > T);
save_coeff = sum(abs(X)>T,2) > 0;
save_coeff = repmat(save_coeff,1,J);
keep_coeff = (thresh_coeff + save_coeff) > 0;

X_hat_indep = X_n .* thresh_coeff; 
X_hat_joint = X_n .* keep_coeff;

risk_indep = norm(X_hat_indep(:,1) - X(:,1))^2
risk_joint = norm(X_hat_joint(:,1) - X(:,1))^2
