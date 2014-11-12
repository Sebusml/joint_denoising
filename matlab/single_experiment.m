% Repeat the result shown in the poster/thesis proposal

clear all
close all
clc

N = 1024;
J = 2:50;
K = 50;
sigma = 0.4;

% N = 4096;
% J = 20;
% K = 50;
% sigma = 0.2;

risks_indep = zeros(length(J), 1);
risks_joint = zeros(length(J), 1);
risks_oracle = zeros(length(J), 1);
improvement_theoretical = zeros(length(J), 1);
improvement_simulation = zeros(length(J), 1);

for i = 1:length(J),
    [~, ~, ~, risk_indep, risk_joint, risk_ub, risk_oracle, improv] = ...
                        joint_denoising_time_2(J(i), N, K, sigma);
    risks_indep(i) = mean(risk_indep);
    risks_joint(i) = mean(risk_joint);
    risks_oracle(i) = mean(risk_oracle);
    improvement_theoretical(i) = mean(improv);
    improvement_simulation(i) = mean(risk_indep) - mean(risk_joint);
end

sigma_x = 1;
T = sigma * sqrt(2*log(N));
improvement_theoretical = improvement(T, J, sigma_x, sigma, K, N);

%% Plotting
do_plot = false;
if do_plot,
    figure
    plot(J, risks_indep, J, risks_joint, J, risks_oracle,'r')
    legend('indep', 'joint', 'oracle')

    figure
    plot(J, improvement_simulation, J, improvement_theoretical)
    legend('simulation', 'theoretical')
end

%% Save data
save_data = true;
if save_data,
    save ('joint_denoising_alt1', '-V7')
end