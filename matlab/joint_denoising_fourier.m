% Joint denoising in the fourier domain

clear all
close all
clc

J = 10;
N = 2048;
K = 5; %round(0.1*N);
X = make_jsm2_signals(J,N,'fourier',K);


% Add noise
sigma = 0.002;
noise = sigma * randn(N,J);
X_n = X + noise;

% Get the Fourier coeff
fc = fft(X);
fc_n = fft(X_n);

% Threshold denoising
fc_hat = fc_n;
T = sigma*sqrt(2*log(N));
for i = 1:N,
    for j = 1:J,
        if abs(fc_n(i,j)) < T,
            fc_hat(i,j) = 0;
        end
    end
end
X_hat = real(ifft(fc_hat));


figure, plot(X(:,1)), hold on, plot(X_n(:,1),'g'), plot(X_hat(:,1),'r')
figure, plot(abs(fft(X_n(:,1))))
figure, plot(abs(fft(X_hat(:,1))))
SNR(X(:,1),X_n(:,1))