clear all
close all
clc

S = 1000;
T = 0.4;
sigma_x = 4;

trials = 1000;
J = 2:1:10;

sums = zeros(trials,1);
N = zeros(trials,1);
sum_means = zeros(length(J),1);
N_mean = zeros(length(J),1);
for j = 1:length(J),
    for n = 1:trials,
        X = abs(sigma_x*randn(S,J(j)));
        res = sum(abs(X(:, 1:J(j)-1)) > T, 2) > 0;
        res =  res .* (abs(X(:,J(j))) < T);
        sums(n) = sum(X(res == 1,J(j)).^2);
        N(n) = sum(res);
        X(res==1,J(j));
    end
    N_mean(j) = mean(N);
    sum_means(j) = mean(sums);
end

% plot(J,sum_means)
% ylabel('Conditional \Sigma \theta_k')
% xlabel('J')
% figure
% plot(J,N_mean)
% xlabel('J')

% sum(res)/S
p1 = normcdf(T, 0, sigma_x) - normcdf(-T, 0, sigma_x);
p = (p1*(1-p1.^(J-1)));

phi = @(x) exp(-0.5*x.^2)/sqrt(2*pi);

E_theta_theo = sigma_x^2*(1 - 2*(T/sigma_x)*phi(T/sigma_x)/p1);

sum_theo = p*S*E_theta_theo;


figure, hold on
plot(J,sum_means), plot(J,sum_theo,'rx')
ylabel('Conditional \Sigma \theta_k')
xlabel('J')
legend('Monte Carlo', 'Theorethical','Location','NorthWest')
% figure
% plot(J,N_mean)
% xlabel('J')



% clc 
% T = 0.4;
% J = 3;
% X = randn(50,J)
% res_1 = sum(abs(X(:, 1:J-1)) > T, 2) > 0;
% res_2 = abs(X(:,J)) < T;
% res_3 = res_1 .* res_2;
% 
% res = sum(abs(X(:, 1:J-1)) > T, 2) > 0;
% res = res .* (abs(X(:,J)) < T);
% 
% [res_1 res_2 res_3 res]
% sum(res)
