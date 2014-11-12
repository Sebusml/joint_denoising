
clear all
close all
clc

N = 10000;
sigma_x = 1;
sigma_y = 1.2;
T = 1;
x = sigma_x * randn(N,1);
y = sigma_y * randn(N,1);

u = x + y;
v = x - y;


z = 0;
for i = 1:N,
     if abs(x(i) + y(i)) < T,
         z(end+1) = x(i)^2 - y(i)^2;
     end
end

z_mean = mean(z)


rho = (sigma_x^2 - sigma_y^2) / (sigma_x^2 + sigma_y^2);
sig_u = sqrt(sigma_x^2 + sigma_y^2);
phi = @(x) exp(-0.5*x.^2)/sqrt(2*pi);
p = normcdf(T, 0, sig_u) - normcdf(-T, 0, sig_u);

z_mean_theo = rho * sig_u^2 * (1 - 2*(T/sig_u)*phi(T/sig_u)/p) 



% zi = x.^2 - y.^2;
% 
% mean(zi)
% mean(z)
% 
% xmax = 15;
% subplot(211)
% hist(zi,20);
% h = gca;
% title('x^2 - y^2')
% subplot(212)
% hist(z,10)
% xlim(get(h,'XLim'))
% title(['x^2 - y^2 | abs(x+y) <'  num2str(T)])



% Mixed chi square
% a1 = 6; a2 = -2;
% z = a1*x.^2 + a2*y.^2;
% m = mean(z)
% m_theo = a1+a2

% Is the pdf of x^2 different than the pdf of X^2| abs(x+y) < T ?
% z1 = x.^2;
% z2 = 0;
% for i = 1:N,
%      if abs(x(i) + y(i)) < T,
%          z2(end+1) = x(i)^2;
%      end
% end
% 
% m_z1 = mean(z1)
% m_z2 = mean(z2)
% 
% chi21 = @(x) x.^-.5.*exp(-x/2)/sqrt(2*pi).*(x>=0);
% subplot(211), hist(z1,15), title('Not conditioned'), h = gca;
% hold on, x = linspace(0,16); plot(x, N*chi21(x),'r')
% subplot(212), hist(z2,15), title('Conditioned'), xlim(get(h,'XLim'))
% hold on, x = linspace(0,16); plot(x, length(z2)*chi21(x),'r')


% chi21 = @(x) x.^-.5.*exp(-x/2)/sqrt(2*pi).*(x>=0);
% N = 10e3;
% Z = randn(N,1).^2;
% hist(Z,15), hold on
% x = linspace(-0.5,20,100); 
% plot(x, N*chi21(x))

% 
% phi = @(x) exp(-0.5*x.^2)/sqrt(2*pi);
% chi21 = @(x) x.^.5.*exp(-x/2)/sqrt(2*pi).*(x>=0);
% close all
% N = 10e3;
% x = randn(N,1);
% y = randn(N,1);
% z = x.^2 + y.^2;
% hist(x,10);
% hold on
% x = linspace(-4,4,100); 
% plot(x, N*phi(x))

% hist(z,20);
% hold on
% x = linspace(0,16,100); 
% plot(x,N*chi22(x),'r')
% 
% 
% f = @(x,y) (x.*y).^-.5.*exp(-(x+y)/2)/2/pi;
