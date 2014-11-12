function imp = improvement(T, J, sigma_x, sigma_w, S, N)

sigma_y = sqrt(sigma_x.^2 + sigma_w.^2); % Measurement variance
p1 = normcdf(T, 0, sigma_y) - normcdf(-T, 0, sigma_y);
p = (p1.*(1-p1.^(J-1)));
p2 = normcdf(T, 0, sigma_y) - normcdf(-T, 0, sigma_y);
rho = (sigma_x.^2 - sigma_w.^2) / (sigma_x.^2 + sigma_w.^2);
phi = @(x) exp(-0.5*x.^2)/sqrt(2*pi);

E_theta_theo = rho .* sigma_y.^2 .* (1 - 2*(T./sigma_y).*phi(T./sigma_y)./p2);

imp = p.*S.*E_theta_theo;


% Add the (negative) improvement due to occurrence of "condition 2"

p2 = normcdf(T/sigma_w) - normcdf(-T/sigma_w);
E_card = (N-S) * p2 * (1 - p2.^(J-1));
E_cond = sigma_w.^2 * (1 - (2*T*phi(T/sigma_w))/(sigma_w*p2));
imp = imp - E_card * E_cond;
