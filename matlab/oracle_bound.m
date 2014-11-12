function ro = oracle_bound(X, sigma)

mins = min(X.^2, repmat(sigma^2, size(X,1), size(X,2)));
ro = sum(mins);