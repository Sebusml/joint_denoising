clear all
clc

N = 1024;
S = 102;
J = 2:30000;
s_w = 0.2;
T = s_w * sqrt(2*log(N));
phi = @(x) exp(-0.5*x.^2)/sqrt(2*pi);
p2 = normcdf(T/s_w) - normcdf(-T/s_w);
E_card = (N-S) * p2 * (1 - p2.^(J-1));
E_cond = s_w^2 * (1 - (2*T*phi(T/s_w))/(s_w*p2));
imp = E_card * E_cond;

plot(J,imp)
