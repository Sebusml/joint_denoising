clear all
clc


% Experiment 1
N = 4096;
sparsity_rate = 0.05:0.05:0.8;
K = round(N*sparsity_rate);
sigma = linspace(0.05,0.5,25);
J = 2:1:50;


% Experiment 2
% N = 4096;
% sparsity_rate = [0.05 0.1];
% K = round(N*sparsity_rate);
% sigma = [0.4 0.5];
% J = 2:5:500;


tic
data = run_joint_denoising_experiment(N, K, sigma, J);
toc


% data(5:6)

save data
