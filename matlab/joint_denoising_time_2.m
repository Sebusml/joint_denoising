function [SNR_original, SNR_indep, SNR_joint, risk_indep, risk_joint, ...
    upper_bound_indep, risk_oracle, improv] = ...
    joint_denoising_time_2(J,N,K,sigma)
% Joint denoising in the time domain

X = make_jsm2_signals(J,N,'time',K,'normal');


% Add noise
noise = sigma * randn(N,J);
X_n = X + noise;

T = sigma * sqrt(2*log(N));


thresh_coeff = (abs(X_n) > T);
save_coeff = sum(abs(X_n)>T,2) > 0;
save_coeff = repmat(save_coeff,1,J);
keep_coeff = (thresh_coeff + save_coeff) > 0;

X_hat_indep = X_n .* thresh_coeff; 
X_hat_joint = X_n .* keep_coeff;


SNR_original = zeros(J,1);
SNR_indep = zeros(J,1);
SNR_joint = zeros(J,1);
risk_indep = zeros(J,1);
risk_joint = zeros(J,1);
improv = zeros(J,1);
% risk_joint_alt = zeros(J,1);
for j = 1:J,
    SNR_original(j) = norm(X(:,j)) / norm(X(:,j) - X_n(:,j));
    SNR_indep(j) = norm(X(:,j)) / norm(X(:,j) - X_hat_indep(:,j));
    SNR_joint(j) = norm(X(:,j)) / norm(X(:,j) - X_hat_joint(:,j));
    risk_indep(j) = (norm(X(:,j) - X_hat_indep(:,j)))^2;
    risk_joint(j) = (norm(X(:,j) - X_hat_joint(:,j)))^2;
    sigma_x = sqrt(var(X(X(:,j)~=0,j))) ; 
    improv(j) = improvement(T, J, sigma_x, sigma, K,N);
    %risk_joint_alt(j) = (norm(X(:,j) - X_hat_joint_alt(:,j)))^2;
end

upper_bound_indep = risk_upper_bound(X, sigma);
risk_oracle = oracle_bound(X, sigma);

%% Save data
save_data = true;
if save_data,
    save ('jd_single_signal', '-V7')
end
%%%%%%%%% Plots %%%%%%%%%%
% ns  = 1:N;
% m = max(abs(X_n(:,1))) + 0.4;
% figure(1)
% xx = X(:,1);
% plot(ns(xx~=0),xx(xx~=0),'bo','MarkerSize',4)
% hold on
% plot(xx)%,'linewidth',2)
% axis([1 N -m m])
% xlabel('$k$')
% title(['Original signal S = ' num2str(K)])
% h1 = line([1 N],[T T], 'color','r');
% h2 = line([1 N],[-T -T], 'color','r');
% % annotate(h1,1000,'$T$','ur',pi/4)
% % annotate(h2,1000,'$-T$','lr',pi/4)
% 
% figure(2)
% plot(X_n(:,1),'linewidth',1)
% axis([1 N -m m])
% xlabel('$k$')
% title(['Noisy signal $\sigma_w$ = ' num2str(sigma)])
% h1 = line([1 N],[T T], 'color','r');
% h2 = line([1 N],[-T -T], 'color','r');
% % annotate(h1,1000,'$T$','ur',pi/4)
% % annotate(h2,1000,'$-T$','lr',pi/4)
% 
% figure(3)
% clf
% set(gcf,'position',[10 10 500 400])
% set(gcf,'papersize',[5 4])
% set(gcf,'paperposition',[0 0 5 4])
% xj = X_hat_joint(:,1);
% xi = X_hat_indep(:,1);
% plot(ns, xj,'b','linewidth',1)
% hold on
% plot(ns, xi,'r','linewidth',1)
% legend('Joint','Independent')
% plot(ns(xj~=0),xj(xj~=0),'bo','MarkerSize',4)
% plot(ns(xi~=0),xi(xi~=0),'r.')
% axis([1 N -m m])
% xlabel('$k$')
% h1 = line([1 N],[T T], 'color','r');
% h2 = line([1 N],[-T -T], 'color','r');
% annotate(h1,1000,'$T$','ur',pi/4)
% annotate(h2,1500,'$-T$','lr',pi/4)
% plotpdftex(1,'sym_1',[1 1])
% plotpdftex(2,'sym_2',[1 1])
% plotpdftex(3,'sym_3',[1 1])

%%%%%%%%%%% End plots %%%%%%%%%%%%%%%

% s = sprintf('Average improvement: %2.2f', mean(SNR_joint) - mean(SNR_indep));
% disp(s)

% Old and slow way

% Independent denoising
% X_hat_indep = X_n;
% for i = 1:N,
%     for j = 1:J,
%         if abs(X_n(i,j)) < T,
%             X_hat_indep(i,j) = 0;
%         end
%     end
% end

% X_hat_indep = X_n .* (abs(X_n) > T);

% Joint denoising
% X_hat_joint = X_n;
% for i = 1:N,
%     for j = 1:J,
%         if abs(X_n(i,j)) < T,
%             % look if at least one of the other coefficients if bigger than
%             % T
%             n_big = sum(abs(X_n(i,:))>T);
%             if n_big == 0,
%                 X_hat_joint(i,j) = 0;
%             end
%         end
%     end
% end

% I am not using the alternative method
% Joint denoising alternative method
% X_hat_joint_alt = X_n;
% for i = 1:N,
%     for j = 1:J,
%         if abs(X_n(i,j)) < T,
%             % look if at least one of the other coefficients if bigger than
%             % T
%             n_big = sum(abs(X_n(i,:)) > T);
%             if n_big == 0,
%                 X_hat_joint_alt(i,j) = 0;
%             elseif abs(X_n(i,j)) < 1.5*sigma,
%                 X_hat_joint_alt(i,j) = 0; % even if we know it is a coefficient
%                                           % if it is too small is better to
%                                           % kill it
%             end
%             
%         end
%     end
% end
