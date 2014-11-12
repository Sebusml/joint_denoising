% Joint denoising

clear all
close all
clc

J = 50;
N = 2048;
X = make_jsm2_signals(J,N,'wavelet',10);

L=5; % Bigger L => higher bandwith of the 
% qmf = MakeONFilter('Symmlet',4);
qmf = MakeONFilter('Daubechies',4);
% Add noise
sigma = 0.004;
noise = sigma * randn(N,J);
X_n = X + noise;

% comp_wavelet_support(X,L,qmf,1e-3);
% comp_wavelet_support(X_n,L,qmf,1e-3);

wc = FWT_PO(X,L,qmf);
wc_n = FWT_PO(X_n,L,qmf);

% for i = 1:J,
%     subplot(J,2,2*i-1)
%     PlotWaveCoeff(wc(:,i),L,0)
%     subplot(J,2,2*i)
%     PlotWaveCoeff(wc_n(:,i),L,0)
% end

wc_hat_joint = wc_n;
T = sigma*sqrt(2*log(N));
for i = 2^L+1:N,
    for j = 1:J,
        if abs(wc_n(i,j)) < T,
            wc_hat_joint(i,j) = 0;
        end
    end
end
X_hat_indep = IWT_PO(wc_hat_joint,L,qmf);

% Joint denoising
wc_hat_joint = wc_n;
for i = 2^L+1:N,
    for j = 1:J,
        if abs(wc_n(i,j)) < T,
            % look if at least one of the other coefficients if bigger than
            % T
            n_big = sum(abs(wc_n(i,:))>T);
            if n_big == 0,
                wc_hat_joint(i,j) = 0;
            end
        end
    end
end
X_hat_joint = IWT_PO(wc_hat_joint,L,qmf);

% Trying my other approach
% wc_mean = mean(wc_n,2);
% X_mean = mean(X_n,2);
% noise_mean = X_mean - mean(X,2);

% Try thresholding the mean of the noisy signals
% wc_mean_hat = wc_mean;
% count = 0;
% for i = 2^L+1:N,
%      if abs(wc_mean_hat(i)) < T/sqrt(J),
%          wc_mean_hat(i) = 0;
%          count = count +1 ;
%      end
% end
% X_mean_hat = IWT_PO(wc_mean_hat,L,qmf);
% SNR(mean(X,2), X_mean_hat)
% SNR(mean(X,2), X_mean)
% figure, hold on, plot(X_mean), plot(X_mean_hat,'r')
% 
% supp_hat = abs(wc_mean_hat) > T/sqrt(J);

% supp = abs(mean(wc,2))>0.1e-3;
% % Test to get the right support
% % wc_test = wc .* repmat(supp,1,J);
% % X_test = IWT_PO(wc_test,L,qmf);
% % norm(X-X_test,'fro')
% % sum(supp)
% % figure,plot(X(:,1)),hold on, plot(X_test(:,1),'r')
% sum(abs(supp-supp_hat))
% 
% wc_hat_2 = wc_n .* repmat(supp_hat, 1, J);
% X_hat_joint_2 = IWT_PO(wc_hat_2,L,qmf);

% figure, hold on
% plot(X(:,1),'g'), plot(X_hat_joint_2(:,1),'r'), plot(X_hat(:,1),'b')

%figure, hold on
% plot(wc(:,1),'g'), plot(wc_hat_2(:,1),'rx'), plot(wc_hat(:,1),'b+')
%plot(wc(2^L+1:end,1),'g'), plot(mean(wc_n(2^L+1:end,:),2),'rx'), plot(wc_n(2^L+1:end,1),'b+')

disp('Singnal to noise ratios')
SNR_original = zeros(J,1);
SNR_indep = zeros(J,1);
SNR_joint = zeros(J,1);
for j = 1:J,
    SNR_original(j) = SNR(X(:,j),X_n(:,j));
    SNR_indep(j) = SNR(X(:,j),X_hat_indep(:,j));
    SNR_joint(j) = SNR(X(:,j),X_hat_joint(:,j));
    s = sprintf('Original: %2.2f Independent: %2.2f Joint: %2.2f',...
         SNR_original(j), SNR_indep(j), SNR_joint(j));
    disp(s)
end
s = sprintf('Average improvement: %2.2f', mean(SNR_joint) - mean(SNR_indep));
disp(s)
figure, hold on
plot(SNR_original), plot(SNR_indep), plot(SNR_joint)

figure,hold on
plot(X(:,1),'g'), plot(X_hat_joint(:,1),'r'), plot(X_hat_indep(:,1),'b')
% legend('noisless','joint','independent')

% subplot(311)
% PlotWaveCoeff(wc(:,1),L,0)
% subplot(312)
% PlotWaveCoeff(wc_hat_2(:,1),L,0)
% subplot(313)
% PlotWaveCoeff(wc_hat(:,1),L,0)
% 
% figure
% scale = 20;
% subplot(311)
% PlotWaveCoeff(wc(:,1),L,scale)
% subplot(312)
% PlotWaveCoeff(wc_n(:,1),L,scale)
% subplot(313)
% PlotWaveCoeff(mean(wc_n,2),L,10)
