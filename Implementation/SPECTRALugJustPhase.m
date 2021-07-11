% Y is in R^{L x N} is N observations of signal on Z_L 
% N is the number of noisy signals to sample
% sigma is the standard deviation of the additive Gaussian noise
% for example:
%   SDPugJustPhase([1;2;3], 4, 0.5)
function z = SPECTRALugJustPhase(Y, sigma, DEBUG)
    % distance function of samples: D(i,j,l) = ||y(i)[-l] - y(j)||^2
    % with Discrete Fourier Transform Dhat
	[L N] = size(Y);
    IP = zeros(N,N,L);
    for i = 1:N
        for j = 1:N     % or use symmetry with j = i+1:N
            for l = 1:L
                IP(i,j,l) = circshift(Y(:,i), -l+1)'*Y(:,j);
            end
        end
    end
    
    % We can relax to a semidefinite program as follows. 
    % Minimizing sum_ij D(i,j) = sum_ij sum_k Dhat(i,j)G(@ik,@jk)/L
    % is similar to the SDP:
    % min_{l_1, ..., l_n in [L]} tr(AG)
   
    index = @(i,k) (i-1)*L + k; % index into 4-D matrix G[i,k1][j,k2]   
    
    A = zeros(N*L, N*L);
    for i = 1:N
        for j = 1:N
            for k = 1:L
                for l = 1:L
                    A(index(i,l),index(j,k)) = IP(i,j,mod(l-k,L)+1);% -(norm(y(:,i))^2 + norm(y(:,j))^2);
                end
            end
        end
    end
    
   
     G = A;
    
    
    
    [U,S,U2] = svds(G,2*L);
    
    % VV^T = G. columns are relaxed variables in SDP, V(:,index(i,k))
    %V = U(:, 1:L+1) * sqrt(S(1:L+1, 1:L+1))
    V = U(:, 1:L) * sqrt(S(1:L, 1:L));
    
    basis_vector = zeros(L,1);
    basis_vector(1) = 1;
    rounding = V * (V(0*N+1:0*N+L, 1:L) \ basis_vector);
    
    % print this rounded vector in nice form
    if DEBUG
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
    
end