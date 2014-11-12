clear all
close all
clc

% load d

load data
d = data;

% load data/data_1
% d = data;
% load data/data_2
% d = [d; data];
% load data/data_3
% d = [d; data];
% load data/data_4
% d = [d; data];
% load data/data_5
% d = [d; data];
% load data/data_6
% d = [d; data];
% load data/data_7
% d = [d; data];
% load data/data_8
% d = [d; data];


get_K = @(N, s, J, d) d(d(:,1)==N & d(:,3)==s & d(:,4)==J, [2 5:end]);
get_sigma = @(N, K, J, d) d(d(:,1)==N & d(:,2)==K & d(:,4)==J, [3 5:end]);
get_J = @(N, K, s, d) d(d(:,1)==N & d(:,2)==K & d(:,3)==s, [4 5:end]);
get_J_sigma = @(N, K, d) d(d(:,1)==N & d(:,2)==K, [3 4 5:end]);

% sparsity_rate = 0.05:0.05:0.8;
K = round(N*sparsity_rate);
%Kg = K(5);
% fprintf('Sparsity ratio: %f', Kg/N)

for i = 1:length(K),
    Kg = K(i);
    risks = get_J_sigma(N,Kg,d);
    ri = risks(:,3);
    rj = risks(:,4);
    rb = risks(:,5);
    ro = risks(:,6);
    Ri = zeros(length(J), length(sigma));
    Rj = zeros(size(Ri));
    Rb = zeros(size(Ri));
    Ro = zeros(size(Ri));
    for s = 1:length(sigma),
        for j = 1:length(J),
            idx = risks(:,1) == sigma(s) & risks(:,2) == J(j);
            Ri(j,s) = ri(idx);
            Rj(j,s) = rj(idx);
            Rb(j,s) = rb(idx);
            Ro(j,s) = ro(idx);
        end
    end

    figure
    set(gcf,'papersize',[50 40])
    set(gcf,'paperposition',[0 0 50 40])
    hold on
    mesh(sigma, J, Ri, 'EdgeColor', 'b')

    mesh(sigma, J, Rj, 'EdgeColor', 'r')
    % mesh(sigma, J, Rb, 'EdgeColor', 'g')
    mesh(sigma, J, Ro, 'EdgeColor', 'c')
    xlabel('$\sigma$')
    ylabel('$J$')
    zlabel('$risk$')
    title(sprintf('Sparsity ratio $%.2f$', Kg/N))
    % legend('Independent', 'Joint', 'Bound','Oracle')
    legend('Independent', 'Joint','Oracle')
    % view(50,20)
    %view(35,35)
    view(40,30)

    file_name = ['figs/risk3d_' num2str(Kg)];
    plotpdftex(gcf, file_name)
%     file_name = ['figs/risk3d_' num2str(Kg) '.png'];
%     exportfig(gcf,file_name,'Format','png','Color','cmyk')
    
%     file_name = ['figs/risk3d_' num2str(Kg) '.eps'];
%     exportfig(gcf,file_name,'Format','eps','Color','cmyk')
    
    close
end

for i = 1:length(K),
    Kg = K(i);
    for j = 1:5:length(sigma),
        sigma_w = sigma(j);
        sigma_x = 1;
        T = sigma_w*sqrt(2*log(N));
        risks = get_J(N,Kg,sigma_w,d);
        Js = risks(:,1);
        ri = risks(:,2);
        rj = risks(:,3);
        ro = risks(:,5);

        imp_theo = risks(:,6);
        %imp_theo = improvement(T, J, sigma_x, sigma_w,  Kg);
        imp_sim = ri - rj;
        
        
        
%         subplot(211)
        figure(1)
        set(gcf,'position',[10 10 600 300])
        set(gcf,'papersize',[6 3])
        set(gcf,'paperposition',[0 0 6 3])
        plot(Js,ri,Js,rj,Js,ro);
        xlabel('$J$')
        ylabel('$risk$')
        legend('independent', 'joint', 'oracle','location','east')%, 'joint alt')
        title(sprintf('Sparsity ratio: $%.2f$ Noise variance: $%.3f$',Kg/N, sigma_w))
%         subplot(212)
        figure(2)
        set(gcf,'position',[10 10 600 300])
        set(gcf,'papersize',[6 3])
        set(gcf,'paperposition',[0 0 6 3])
        set(gcf,'papersize',[6 3])
        set(gcf,'paperposition',[0 0 6 3])
        plot(Js, imp_sim, Js, imp_theo);
        xlabel('$J$')
        ylabel('$risk \quad improvement$')
        legend('simulation', 'theoretical','location','southeast')

        file_name = sprintf('figs/risk_%d_%d', Kg, int16(sigma_w*1000));
        plotpdftex(1, [file_name '_a'],[1 1])
        plotpdftex(2, [file_name '_b'], [1 1])
%          file_name = sprintf('figs/risk_%d_%d.png', Kg, int16(sigma_w*1000));
%          exportfig(gcf,file_name,'Format','png','Color','cmyk','Bounds','loose')
%         file_name = sprintf('figs/risk_%d_%d.eps', Kg, int16(sigma_w*1000));
%         exportfig(gcf,file_name,'Format','eps','Color','cmyk');%,'Bounds','loose')
        close

    end
end
