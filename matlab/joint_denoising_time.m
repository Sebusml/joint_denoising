% Joint denoising in the time domain
% Toy example

clear all
close all
clc

J = 3;
N = 1024;
X = zeros(N,J);

n1 = 250;
n2 = 650;
X(n1,1) = 10;
X(n2,1) = 18;
X(n1,2) = 15;
X(n2,2) = 2;
X(n1,3) = 12;
X(n2,3) = 8;

% Add noise
sigma = 1;
noise = sigma * randn(N,J);
X_n = X + noise;

T = sigma * sqrt(2*log(N));

% Independent denoising
X_hat_indep = X_n;
for i = 1:N,
    for j = 1:J,
        if abs(X_n(i,j)) < T,
            X_hat_indep(i,j) = 0;
        end
    end
end

% Joint denoising
X_hat_joint = X_n;
for i = 1:N,
    for j = 1:J,
        if abs(X_n(i,j)) < T,
            % look if at least one of the other coefficients if bigger than
            % T
            n_big = sum(abs(X_n(i,:))>T);
            if n_big == 0,
                X_hat_joint(i,j) = 0;
            end
        end
    end
end

subplot(2*J,2,1)
for j = 1:J,
subplot(2*J,2,2*j-1)
plot(X(:,j))
axis([0 N-1 -2 20])
line([0 N-1], [T T],'color','r','LineStyle',':')

subplot(2*J,2,2*j)
plot(X_n(:,j))
axis([0 N-1 -2 20])
line([0 N-1], [T T],'color','r','LineStyle',':')

subplot(2*J,2,2*j+5)
plot(X_hat_indep(:,j))    
axis([0 N-1 -2 20])
line([0 N-1], [T T],'color','r','LineStyle',':')

subplot(2*J,2,2*j+6)
plot(X_hat_joint(:,j))    
axis([0 N-1 -2 20])
line([0 N-1], [T T],'color','r','LineStyle',':')
end

