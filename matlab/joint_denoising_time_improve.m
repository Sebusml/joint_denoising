function [risk_indep, risk_joint,  risk_joint_improv, risk_oracle] = ...
    joint_denoising_time_improve(J,N,K,sigma)
% Joint denoising in the time domain

X = make_jsm2_signals(J,N,'time',K,'normal');


% Add noise
noise = sigma * randn(N,J);
X_n = X + noise;

T = sigma * sqrt(2*log(N));


thresh_coeff = (abs(X_n) > T);
save_coeff = sum(abs(X_n)>T,2) > 0;
save_coeff = repmat(save_coeff,1,J);
keep_coeff = (thresh_coeff + save_coeff) > 0;

if J < 20,
    Nj = 1;
else
    Nj = 2;
end

%Nj = ceil(sigma * sqrt(2*log(J)));

save_coeff_improv = sum(abs(X_n)>T,2) >= Nj;
save_coeff_improv = repmat(save_coeff_improv,1,J);
keep_coeff_improv = (thresh_coeff + save_coeff_improv) > 0;

X_hat_indep = X_n .* thresh_coeff; 
X_hat_joint = X_n .* keep_coeff;
X_hat_joint_improv = X_n .* keep_coeff_improv;

risk_indep = zeros(J,1);
risk_joint = zeros(J,1);
risk_joint_improv = zeros(J,1);

for j = 1:J,
    risk_indep(j) = (norm(X(:,j) - X_hat_indep(:,j)))^2;
    risk_joint(j) = (norm(X(:,j) - X_hat_joint(:,j)))^2;
    risk_joint_improv(j) = (norm(X(:,j) - X_hat_joint_improv(:,j)))^2;
end

risk_oracle = oracle_bound(X, sigma);