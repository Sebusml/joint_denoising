
clear all
close all
clc

N = 4096;
sparsity_rate = 0.01;
K = round(N*sparsity_rate);
sigma = 0.2;
J = 10;

[~, ~, ~, risk_indep, risk_joint, risk_ub, risk_oracle, improv] = ...
      joint_denoising_time_2(J, N, K, sigma);
