clear all
close all
clc


% Monte Carlo simulation to test the conditional probability of a  random 
% variable. The random variable is normal, and the condition is that its
% absolute value is smaller than T.

N = 2^14;
sigma_w = 0.02;
sigma_x = 1;
T = sigma_w * sqrt(2*log(N));

% Unconditional probability
X = sigma_x * randn(N,1);
p_hat = sum(abs(X)<T)/N;
p_theo = normcdf(T, 0, sigma_x) - normcdf(-T, 0, sigma_x)


% Conditional probability

% Generate a random vector with N entries such that abs(x) < T
X = zeros(N,1);
n = 0;
while 1,
    tmp = sigma_x * randn(N,1); 
    tmp = tmp(abs(tmp) < T);
    n_new = length(tmp);
    if n + n_new >= N,
        X(n+1:N) = tmp(1:N-n);
        break
    else
        X(n+1:n+n_new) = tmp;
        n = n + n_new;
    end
end

E_hat = mean(X.^2)

phi = @(x) exp(-0.5*x.^2)/sqrt(2*pi);

p1 = p_theo;
E_theo = sigma_x^2*(1 - 2*(T/sigma_x)*phi(T/sigma_x)/p1)
% E_theo = sigma_x^2*(1 + 2 * inv(N)^(sigma_w^2/sigma_x^2) / sqrt(2*pi) / p1)
