% Joint denoising in the time domain

clear all
close all
clc

J = 50;
N = 4096;
K = round(0.05*N);
X = make_jsm2_signals(J,N,'time',K,'normal');

% Add noise
sigma = 0.237;
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

% disp('Signal to noise ratios')
SNR_original = zeros(J,1);
SNR_indep = zeros(J,1);
SNR_joint = zeros(J,1);
risk_indep = zeros(J,1);
risk_joint = zeros(J,1);
for j = 1:J,
%     SNR_original(j) = SNR(X(:,j),X_n(:,j));
%     SNR_indep(j) = SNR(X(:,j),X_hat_indep(:,j));
%     SNR_joint(j) = SNR(X(:,j),X_hat_joint(:,j));
    risk_indep(j) = norm(X(:,j) - X_hat_indep(:,j));
    risk_joint(j) = norm(X(:,j) - X_hat_joint(:,j));
    s = sprintf('Original: %2.2f Independent: %2.2f Joint: %2.2f',...
         SNR_original(j), SNR_indep(j), SNR_joint(j));
     disp(s)
end
% s = sprintf('Average improvement: %2.2f', mean(SNR_joint) - mean(SNR_indep));
% disp(s)
% figure, hold on
% plot(SNR_original,'g'), plot(SNR_indep,'b'), plot(SNR_joint,'r')
% legend('original','independent','joint')


figure, hold on
plot(X(:,1),'gx'), plot(X_hat_joint(:,1),'r'), plot(X_hat_indep(:,1),'b')
legend('noisless','joint','independent')

risk_bound = risk_upper_bound(X, sigma);
figure, hold on
plot(risk_bound,'g'), plot(risk_indep,'b'), plot(risk_joint,'r')
legend('upper bound','independent','joint')