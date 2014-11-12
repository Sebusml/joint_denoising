clear all
close all
clc

load data
d = data;

get_K = @(N, s, J, d) d(d(:,1)==N & d(:,3)==s & d(:,4)==J, [2 5:end]);
get_sigma = @(N, K, J, d) d(d(:,1)==N & d(:,2)==K & d(:,4)==J, [3 5:end]);
get_J = @(N, K, s, d) d(d(:,1)==N & d(:,2)==K & d(:,3)==s, [4 5:end]);
get_J_sigma = @(N, K, d) d(d(:,1)==N & d(:,2)==K, [3 4 5:end]);

K = 205;

Kg = 205;
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
mesh(sigma, J, Ro, 'EdgeColor', 'c')
xlabel('$\sigma$')
ylabel('$J$')
zlabel('$risk$')
title(sprintf('Sparsity ratio $%.2f$', Kg/N))
legend('Independent', 'Joint','Oracle')
view(40,30)

save_data = true;
if save_data,
    save ('plot3d_singleK', '-V7')
end
