% %% Leer datos
% data = csvread('patient001_s0010.csv',2,0);
% N = 512*3;
% L = 15;
% offset = 0;
% % Tomar ventana de datos
% X = data(1+offset:offset+N,2:L+1);
% for i=1:15
%     figure (i);
%     plot(X(:,i))
%     ylim([-0.8 0.8]);
%     xlim([0 N]);
% end
% %% Wavelet Descomposition un nivel
% %Perform a single-level decomposition of the signal using the db1 wavelet.
% l_s = N;
% [cA1,cD1] = dwt(X(:,1),'db1');
% A1 = idwt(cA1,[],'db1',l_s);
% D1 = idwt([],cD1,'db1',l_s);
% subplot(1,3,1); plot(A1); title('Approximation A1')
% subplot(1,3,2); plot(D1); title('Detail D1')
% subplot(1,3,3); plot(X(:,1)); title('Original')
% 
% %% Wavelet Descomposition
% [C,L] = wavedec(X(:,1),3,'db1');
% [cD1,cD2,cD3] = detcoef(C,L,[1,2,3]);
% % Aproximation level 3
% A3 = wrcoef('a',C,L,'db1',3);
% % To reconstruct the details at levels 1, 2, and 3, from C
% D1 = wrcoef('d',C,L,'db1',1); 
% D2 = wrcoef('d',C,L,'db1',2); 
% D3 = wrcoef('d',C,L,'db1',3);
% 
% figure(1)
% %plots
% subplot(2,3,1); plot(A3);  
% title('Approximation A3') 
% subplot(2,3,2); plot(D1);  
% title('Detail D1') 
% subplot(2,3,3); plot(D2);  
% title('Detail D2') 
% subplot(2,3,4); plot(D3);  
% title('Detail D3')
% 
% %Reconstruct signal
% A0 = waverec(C,L,'db1');
% subplot(2,3,5); plot(A0);  
% title('Reconstructed Signal')
% figure(2)
% % Compare signal
% subplot(2,1,1);plot(X(:,1));title('Original'); axis off 
% subplot(2,1,2);plot(A3);title('Level 3 Approximation'); 
% axis off
% %% Wavelet Descomposition 2
% [C,L] = wavedec(X(:,1),9,'db1');
% [cD1,cD2,cD3,cD4,cD5,cD6,cD7,cD8,cD9] = detcoef(C,L,[1,2,3,4,5,6,7,8,9]);
% % A3 = wrcoef('a',C,L,'db1',3);
% 
% % figure(1)
% % %plots
% % for i=1:9
% %     subplot(3,3,i); plot(wrcoef('d',C,L,'db1',i));   
% % end
% figure(2)
% plot(C)
% 
% 
% 
% %% Wavelet Package
% fs = 1000;
% T = wpdec(X(:,1),9,'db1','shannon');
% [S,T,F] = wpspectrum(T,fs,'plot');
% 
% %% Comparar Leads
% data = csvread('patient001_s0010.csv',2,0);
% N = 512*3;
% L = 15;
% offset = 0;
% % Tomar ventana de datos
% X = data(1+offset:offset+N,2:L+1);
% C = zeros(size(X));
% %Daubachies
% for j=1:1
%     figure (j)
%     for i=1:15
%         [C(:,i),L] = wavedec(X(:,i),9,sprintf('db%d',j));
%         subplot(4,4,i);plot(C(:,i)); xlim([0 N]);ylim([-1 1]);
%         title(sprintf('Lead %d Db%d',i,j))
%         text(300,-0.5,'C=[cA9 cD9 ... cD1]');
%     end
% end
% %% Surf
% surf(C);
% shading interp
% axis tight
% % title('Original Data: 35 days of electrical consumption');
% xlabel('Lead');
% ylabel('Coef');
% ax = gca;
% ax.View = [-13.5 48];

%% Thresholding
patient = load('s0001_rem.mat');
data = patient.val';
N = 512*3;
L = 15;
offset = 0;

% Tomar ventana de datos
X = data(1+offset:offset+N,:);
C = zeros(size(X));
M = C;
Xa = X;
L = zeros(11,15);
P = zeros(15,1);
error = zeros(N,1);
%Iterar en Leads.
for i=1:15
    %Calculo DWT para Lead i
    [C(:,i),L(:,i)] = wavedec(X(:,i),9,'db1');
    C_sort = sort( abs(C(:,i)) ,'descend');
    %Iterar para distintos M-term
    for n=1:200
        M(:,i) = C(:,i);
        %Thresholding
        M(abs(C(:,i)) < C_sort(n) ,i) = 0;
        %IDWT
        Xa(:,i) = waverec( M(:,i),L(:,i),'db1');
        %Calculo error
        error(n) = norm(X(:,i)-Xa(:,i))*norm(X(:,i)-Xa(:,i));
    end
    %plot error n-esimo
    subplot (4,4,i)   
    title(sprintf('Lead %d ',i))
    plot(error)
    xlim([1 200]); %ylim([0 5]);
end
% Sup M* = 50;
%% EXP2
S = zeros(size(M));
for i=1:15
    %Calculo DWT para Lead i
    [C(:,i),L(:,i)] = wavedec(X(:,i),9,'db1');
    C_sort = sort( abs(C(:,i)) ,'descend');
    M(:,i) = C(:,i);
    %Thresholding
    M(abs(C(:,i)) < C_sort(512/4) ,i) = 0;
    M(find(M(:,i)) ,i) = 1;
end

[r,c] = size(M');                           %# Get the matrix size
imagesc((1:c)+0.5,(1:r)+0.5,M');            %# Plot the image
colormap(gray);

%%
%addpath('.')
load('s0010_rem.mat');


