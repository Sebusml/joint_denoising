% Repeat the result shown in the poster/thesis proposal

clear all
close all
clc

N = 1024;
J = 2:1000;
K = 50;
sigma = 0.4;

risks_indep = zeros(length(J), 1);
risks_joint = zeros(length(J), 1);
risks_oracle = zeros(length(J), 1);
risks_joint_improve = zeros(length(J), 1);
improvement_simulation = zeros(length(J), 1);
improvement_simulation_2 = zeros(length(J), 1);

for i = 1:length(J),
    [risk_indep, risk_joint, risk_joint_improve, risk_oracle] = ...
                        joint_denoising_time_improve(J(i), N, K, sigma);
    risks_indep(i) = mean(risk_indep);
    risks_joint(i) = mean(risk_joint);
    risks_joint_improve(i) = mean(risk_joint_improve);
    risks_oracle(i) = mean(risk_oracle);
    improvement_simulation(i) = mean(risk_indep) - mean(risk_joint);
    improvement_simulation_2(i) = mean(risk_indep) - mean(risk_joint_improve);
    if mod(i,100) == 0,
        i
    end
end

sigma_x = 1;
T = sigma * sqrt(2*log(N));
improvement_theoretical = improvement(T, J, sigma_x, sigma, K, N);

%% Plotting
figure
plot(J, risks_indep, J, risks_joint, J, risks_joint_improve, ...
     J, risks_oracle)
legend('indep', 'joint', 'joint_2', 'oracle', 'location', 'best')

figure
plot(J, improvement_simulation, J, improvement_theoretical, ...
     J, improvement_simulation_2)
legend('simulation', 'theoretical', 'simulation_2', 'location', 'best')

%% Save data
save_data = false;
if save_data,
    save ('joint_denoising_alt2', '-V7')
end
