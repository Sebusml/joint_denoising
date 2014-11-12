%% Leer datos
data = csvread('patient001_s0010.csv',2,0);
N = 512*3;
L = 15;
offset = 0;
% Tomar ventana de datos
X = data(1+offset:offset+N,2:L+1);
for i=1:15
    figure (i);
    plot(X(:,i))
    ylim([-0.8 0.8]);
    xlim([0 N]);
end

