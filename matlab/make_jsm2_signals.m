function  X = make_jsm2_signals(J, N, domain, K, varargin)
%MAKE_JSM_SINGNALS - Make n signals satisfying the JSM-2 model
% 
% Syntax:  X = make_jsm2_signals(n, K, domain)
%
% Inputs:
%    J - Number of signals
%    N - length of the signals
%    domain - In wich domain the signal is sparse ('Fourier', 'Time',
%             'Wavelet')
%    K - Sparsity level of the signals
%    varargin - Extra parameters
%               For time domain signals:
%               'normal': The coefficients follow a normal distribution.
%               'uniform': The coefficients follow a uniform distribution
%                          (default).
% Outputs:
%    X - N by J matrix with the signals. 
%
% Example: 
%
% Other m-files required: Wavelab toolbox

% Author: Alejandro Weinstein
% Colorado School of Mines
% email: alejandro.weinstein@gmail.com
% January 2010; Last revision: 2010-01-07

switch lower(domain)
    case 'fourier'
        X = fourier_signal(J, N, K);
    case 'time'
        X = time_signal(J, N, K, varargin);
    case 'wavelet'
        X = wavelet_signal(J, N, K);
    otherwise
        disp('Unsupported domain')
        X = zeros(N,J);
end

function X = fourier_signal(J, N, K)
    X = zeros(N*J, 1);
    t = (0:N-1)';
    %periods = 2.^(9:-1:(9-K+1));
    periods = 2.^(log2(N):-1:log2(N)-K);
    for k = 1:J,
        amp = rand(K,1); % amplitudes
        freq = 2*pi./periods; % frequencies
        phases = 2*pi*rand(K,1);
        for i = 1:K,
            X((1:N) + (k-1)*N) = X((1:N) + (k-1)*N) + ...
                amp(i)*cos(freq(i)*t + phases(i));
        end
        X((1:N) + (k-1)*N) = X((1:N) + (k-1)*N) / norm(X((1:N) + (k-1)*N));    
    end
    X = reshape(X, N, J);
   
function X = time_signal(J, N, K, varargin)
    tmp = randperm(N);
    support = sort(tmp(1:K));
    X = zeros(N,J);
    if strcmp(varargin{1},'uniform'),
        a = 1;
        coeff = 2*a*rand(K,J) - a ;
    else
        mu = 0; % Mean value of the coefficients
        sigma = 1; % Sigma of the coefficients
        coeff = sigma * randn(K,J) + mu;
    end
%     bias = sum(coeff,1) / K;
%     coeff = coeff - repmat(bias,K,1); % Force zero mean
%     coeff_norm = sqrt(diag(coeff' * coeff))';
%     coeff = coeff ./ repmat(coeff_norm,K,1); % Force norm to be 1
    X(support,:) = coeff;
    
function X = wavelet_signal(J, N, K)
    X = zeros(N,J);
    t = ((1:fix(N/5)) ./fix(N/5))';
    sigma = 5;
    for i = 1:J,
        par = sigma * randn(1,12) + [20 4 10 45 40 100 16 8 16 20 4 20]; 
        sig1=par(1)*(t.^3+t.^2+par(2));
        sig2=par(3)*t.^3 + par(4);
        sig3=par(5)*(2.*t.^3+t) + par(6);
        sig4=par(7)*t.^2+par(8)*t+par(9);
        sig5=par(10)*(t+par(11));
        sig6(1:fix(N/10))=ones(1,fix(N/10));
        sig6=sig6*par(12);

        sig(1:fix(N/5))=sig1;
        sig(2*fix(N/5):-1:(fix(N/5)+1))=sig2;
        sig((2*fix(N/5)+1):3*fix(N/5))=sig3;
        sig((3*fix(N/5)+1):4*fix(N/5))=sig4;
        sig((4*fix(N/5)+1):5*fix(N/5))=sig5(fix(N/5):-1:1);
        diff=N-5*fix(N/5);
        sig(5*fix(N/5)+1:N)=sig(diff:-1:1);
        sig((fix(N/20)+1):(fix(N/20)+fix(N/10)))=ones(1,fix(N/10))*10;
        sig((N-fix(N/10)+1):(N+fix(N/20)-fix(N/10)))=ones(1,fix(N/20))*150;

        % zero-meaN
        bias=sum(sig)/N;
        sig = sig - bias;
        sig = sig / norm(sig); % Normalize to norm 1
        X(:,i) = sig;
    end
