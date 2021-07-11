% I have not tested this because I do not have access to cvx.

% Y is in R^{L x N} is N observations of signal on Z_L 
% N is the number of noisy signals to sample
% sigma is the standard deviation of the additive Gaussian noise
% for example:
%   SDPugJustPhase([1;2;3], 4, 0.5)
function z = UGJustPhase_opt2(Y, sigma, DEBUG)
    [L N] = size(Y);
    index  = @(i,k) (i-1)*L + k;  % index into 4-D matrix G[i,k1][j,k2]   
    index2 = @(i,k) (k-1)*N + i; % index into 4-D matrix G[i,k1][j,k2] with N-dimension first   
    
    tic
    
    P = zeros(N*L,N*L);         % permutation matrix
    for i = 1:N
        for k = 1:L
            P(index(i,k),index2(i,k)) = 1;
        end
    end
    DFT = kron(conj(dftmtx(L))/sqrt(L),eye(N))*P';

    %% componentwise fft is much faster for large L
    V_FFT = zeros(N*L,L);   
    % V_FFT is LDL decomposition of coefficient matrix C, so that after 
    % normalization, C = V_FFT*V_FFT'/L.
    for i = 1:N
        Fc = fft(arrayfun(@(k) circshift(Y(:,i), k-1)' * Y(:,1), 1:L));
        for k = 1:L
            V_FFT(index2(i,k),k) = Fc(k)/abs(Fc(k));
        end
    end
    
    %% rounding
    V = DFT' * V_FFT;
    [rL,rN] = deal(1,1);
    basis_vector = zeros(L,1);
    basis_vector(rL) = 1;
    rounding = real(V * (V(index(rN,1):index(rN,L), 1:L) \ basis_vector));
    
    % print this rounded vector in nice form
    if DEBUG == true
        fprintf('\nRounding from the vectors: ');
        for i=1:N
            disp(rounding((i-1)*L+1:(i-1)*L+L)');
        end
    end
    
    Restimated = zeros(N,1);
    z = zeros(L,1); 
    % for each signal, choose rotation with greatest indicator value
    for i = 1:N
        [M, Restimated(i)] = max(rounding(index(i,1):index(i,L)));
        z = z + circshift(Y(:,i), -Restimated(i)+1);
    end
    z = real(PowerSpectrumCorrection(z/N,Y,sigma));   % project result onto real line
    
    fprintf('Total time: %f\n',toc);
end