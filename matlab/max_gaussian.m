% Numerical test of the expected value of the mamimum value of a vector
% of length N, with i.i.d. entries draw from a normal distribution
% Asymptotically, the maximum value is sigma*sqrt(2*ln(N))

clc
I = 50;
N = round(logspace(2,7,I));
sigma = 1;
K = 500;
max_mean = zeros(I,1);
max_hat = zeros(I,1);
for i = 1:length(N),
    max_hat(i) = sigma*sqrt(2*log(N(i)));
    maxs = zeros(K,1);
    for k = 1:K,
        maxs(k) = max(sigma*abs(randn(N(i),1)));
    end
    max_mean(i) = mean(maxs);
end
figure, hold on
semilogx(N, max_hat)
semilogx(N, max_mean,'r')
legend('Asymptotic max', 'Real max',0)
xlabel('N')
ylabel('E(max(w))')
title('Maximum value of vector w with length N; w ~ N(0,1) i.i.d.')