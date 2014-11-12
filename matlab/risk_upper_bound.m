function bound = risk_upper_bound(X, sigma)

n = size(X,1);
ob = oracle_bound(X, sigma);
bound = (2*log(n) + 1) * (sigma^2 + ob);