clear all
clc


% Monte Carlo simulation to test the probability that one random variable 
% is smaller than T, and at least one of n-1 random variables are greater
% than T.
% The random variables are assumed to be i.i.d. with pdf N(0,1).

T = 0.8;
p1 = normcdf(T) - normcdf(-T);
p2 = 1 - p1;
    
trials = 10000;
n = 10;
X = abs(randn(trials,n));

res = zeros(trials,1);
for i = 1:n-1,
    res = res | X(:,i)>T;
end
res = res & X(:,n)<T;

% p_hat = sum(res) / trials % Estimated probability
% 
% p_theo = p1 * (1 - p1^(n-1)) % Theoretic probability


% Monte Carlo simulation to test the number of times that one random variable 
% is smaller than T, and at least one of n-1 random variables are greater
% than T, when trying S times 
% The random variables are assumed to be i.i.d. with pdf N(0,1).

S = 5;
trials = 5000;
N_t = zeros(trials,1);
for j = 1:trials,
    N = zeros(S,1);
    for k = 1:S,
        X = abs(randn(n,1));
        res = 0;
        for i = 1:n-1,
            res = res | X(i) > T;
        end
        res = res & X(n) < T;
        N(k) = res;
    end
    N_t(j) = sum(N);
end

E_N_hat = mean(N_t)
E_N_theo = S*p1*(1-p1^(n-1))
