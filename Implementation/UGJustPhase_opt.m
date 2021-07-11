% I have not tested this because I do not have access to cvx.

% Y is in R^{L x N} is N observations of signal on Z_L 
% N is the number of noisy signals to sample
% sigma is the standard deviation of the additive Gaussian noise
% for example:
%   SDPugJustPhase([1;2;3], 4, 0.5)
function z = UGJustPhase_opt(Y, sigma, DEBUG)
    % distance function of samples: D(i,j,l) = ||y(i)[-l] - y(j)||^2
    % with Discrete Fourier Transform Dhat
	[L N] = size(Y);
    D = zeros(N,N,L);
    index = @(i,k) (i-1)*L + k;  % index into 4-D matrix G[i,k1][j,k2]   
    index2 = @(i,k) (k-1)*N + i; % index into 4-D matrix G2[i,k1][j,k2] with N-dimension first   
    
    %% slow initialization method
    tic
    for i = 1:N
        for j = i:N     % or use symmetry with j = i+1:N
            for l = 1:L
                D(i,j,l) = circshift(Y(:,i), -l+1)' * Y(:,j); % norm(circshift(Y(:,i), -l+1) - Y(:,j))^2;
            end
            for l = 1:L
                D(j,i,l) = D(i,j,mod(1-l,L)+1);
            end
        end
    end
    fprintf('Timing 0: %f\n',toc);
    
    % We can relax to a semidefinite program as follows. 
    % Minimizing sum_ij D(i,j) = sum_ij sum_k Dhat(i,j)G(@ik,@jk)/L
    % is similar to the SDP:
    % min_{l_1, ..., l_n in [L]} tr(AG)
   
    tic
    A = zeros(N*L, N*L);
    for i = 1:N
        for j = 1:N 
            for k = 1:L
                for l = 1:L
                    A(index(i,l),index(j,k)) = D(i,j,mod(l-k,L)+1); % -(norm(y(:,i))^2 + norm(y(:,j))^2)
                end
            end
        end
    end
    fprintf('Timing 1: %f\n',toc);
    
    tic
    
    P = zeros(N*L,N*L); % permutation matrix
    for j = 1:N*L
        i = mod(j-1,N)*L + ceil(j/N); 
        P(i,j) = 1;
    end
    DFT = kron(conj(dftmtx(L))/sqrt(L),eye(N))*P';
    D1 = DFT * A * DFT';
    
    for i = 1:N
        for j = 1:N
            for k = 1:L
                D1(index2(i,k),index2(j,k)) = D1(index2(i,k),index2(j,k)) / abs(D1(index2(i,k),index2(j,k))) / L;
            end
        end
    end
    fprintf('Timing 2: %f\n',toc);
    
    
    
    %% componentwise fft is much faster for large L
    tic
    V_FFT = zeros(N*L,L);   
    % V_FFT is LDL decomposition of coefficient matrix C, so that after 
    % normalization, C = V_FFT*V_FFT'/L.
    for i = 1:N
        Fc = fft(arrayfun(@(k) circshift(Y(:,i), k-1)' * Y(:,1), 1:L));
        for k = 1:L
            V_FFT(index2(i,k),k) = Fc(k)/abs(Fc(k));
        end
        
%   To complete the coefficient matrix C, we may alternatively compute:
%       for j = i:N
%           Fc = fft(arrayfun(@(k) circshift(Y(:,i), k-1)' * Y(:,j), 1:L));
%           for k = 1:L
%               C(i,j,k) = Fc(k)/abs(Fc(k));
%               C(j,i,k) = conj(C(i,j,k));
%           end
%       end
    end
    
    fprintf('Timing 1.5: %f\n',toc);
    
    %% obtain SDP matrix G
    tic
    
    G = DFT' * D1 * DFT;
    G = real(G);    % cleanup precision error on imaginary components
    
    %% rounding
    [U,S,U2] = svds(G,L);
    % VV^T = G. columns are relaxed variables in SDP, V(:,index(i,k))
    V = U(:, 1:L) * sqrt(S(1:L, 1:L));   
    
    basis_vector = zeros(L,1);
    basis_vector(1) = 1;
    rounding = V * (V(0*N+1:0*N+L, 1:L) \ basis_vector);
    
    fprintf('Timing 3: %f\n',toc);
    
    tic
    
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
        best = -1000;
        for k = 1:L
            prod = rounding(index(i,k));
            if prod > best
                Restimated(i) = k - 1;
                best = prod;
            end
        end
        z = z + circshift(Y(:,i), -Restimated(i));
    end
    z = z / N;

    z = PowerSpectrumCorrection(z,Y,sigma);
    
    fprintf('Timing 4: %f\n',toc);
    
    %% rounding 2
    V = DFT' * V_FFT;
    basis_vector = zeros(L,1);
    basis_vector(1) = 1;
    rounding = real(V * (V(0*N+1:0*N+L, 1:L) \ basis_vector));
    
    % print this rounded vector in nice form
    if DEBUG == true
        fprintf('\nRounding from the vectors: ');
        for i=1:N
            disp(rounding((i-1)*L+1:(i-1)*L+L)');
        end
    end
    
end