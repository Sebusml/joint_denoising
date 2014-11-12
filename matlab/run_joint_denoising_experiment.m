function data = run_joint_denoising_experiment(N, K, sigma, J)

len_N = length(N);
len_K = length(K);
len_sigma = length(sigma);
len_J = length(J);
n_rows = len_N * len_K * len_sigma * len_J;
data = zeros(n_rows, 9);

i = 1;
for i_N = 1:len_N,
    for i_K = 1:len_K,
        for i_sigma = 1:len_sigma,
            for i_J = 1:len_J,
                [~, ~, ~, risk_indep, risk_joint, risk_ub, risk_oracle, improv] = ...
                    joint_denoising_time_2(J(i_J), N(i_N), K(i_K), ...
                                           sigma(i_sigma));
%                     data(i, :) = [N(i_N) K(i_K) sigma(i_sigma) J(i_J)...
%                         mean(risk_indep) mean(risk_joint) mean(risk_ub)...
%                         mean(risk_oracle), mean(risk_joint_alt)];
                    data(i, :) = [N(i_N) K(i_K) sigma(i_sigma) J(i_J)...
                        mean(risk_indep) (risk_joint(1)) mean(risk_ub)...
                        mean(risk_oracle), mean(improv)];
                   if mod(i,10) == 0, 
                        i
                   end
                   i = i + 1;
            end
        end
    end
end
