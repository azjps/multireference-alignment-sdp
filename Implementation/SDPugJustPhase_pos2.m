% I have not tested this because I do not have access to cvx.

% Y is in R^{L x N} is N observations of signal on Z_L 
% N is the number of noisy signals to sample
% sigma is the standard deviation of the additive Gaussian noise
% for example:
%   SDPugJustPhase([1;2;3], 4, 0.5)
function z = SDPugJustPhase_pos2(Y, sigma, DEBUG)
    % distance function of samples: D(i,j,l) = ||y(i)[-l] - y(j)||^2
    % with Discrete Fourier Transform Dhat
	[L N] = size(Y);
    D = zeros(N,N,L);
    for i = 1:N
        for j = 1:N     % or use symmetry with j = i+1:N
            for l = 1:L
                D(i,j,l) = circshift(Y(:,i), -l+1)' * Y(:,j); % norm(circshift(Y(:,i), -l+1) - Y(:,j))^2;
            end
        end
    end
    
    % We can relax to a semidefinite program as follows. 
    % Minimizing sum_ij D(i,j) = sum_ij sum_k Dhat(i,j)G(@ik,@jk)/L
    % is similar to the SDP:
    % min_{l_1, ..., l_n in [L]} tr(AG)
   
    index = @(i,k) (i-1)*L + k;  % index into 4-D matrix G[i,k1][j,k2]   
    index2 = @(i,k) (k-1)*N + i; % index into 4-D matrix G2[i,k1][j,k2] with N-dimension first   
    
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
    
    % D
    % A
    
    A;
    % http://www.math.zju.edu.cn/amjcu/B/200401/040103.pdf
    % A is unitarily similar to block diagonal matrix
    % kron(eye(N),dftmtx(L)) * A * kron(eye(N),dftmtx(L))'
    
    P = zeros(N*L,N*L);
    for j = 1:N*L
        i = mod(j-1,N)*L + ceil(j/N); 
        P(i,j) = 1;
    end
    % P
    % P'*A*P
    DFT = kron(dftmtx(L)/sqrt(L),eye(N))*P';
    D1 = DFT * A * DFT';

    cvx_begin sdp
        variable G(N*L, N*L) symmetric
        dual variable Q
        maximize( trace(G*A) )
        % variable x(N,L)
        
        % constraint (1): G[(i,k1)][(i,k2)] = G[(i,k1-k2)][(i,0)]
        for i = 1:N
            for k = 1:L
                G(index(i,k),index(i,k)) == 1/L
            end
            % sum(diag(G(index(i,1):index(i,L), index(i,1):index(i,L)))) == 1
        end
        % constraint (2)
        for i = 1:N
            for l = 1:L
                for k = l+1:L
                    G(index(i,l),index(i,k)) == 0
                end
            end
        end
        % constraint (3): check matrix entries are non-negative
         for i = 1:N
             for j = i+1:N
                 for l = 1:L
                     for k = 1:L
                         0 <= G(index(i,l),index(j,k)) % <= G(index(i,l),index(i,l))
                     end
                 end
             end
         end
        % constraint (4)
        % sum G_{ij} = 1
        for i = 1:N
            for j = i+1:N
                sum(sum(G(index(i,1):index(i,L), index(j,1):index(j,L)))) == 1
            end
        end
        
        G >= 0 : Q;
    cvx_end
    
    for i = 1:N*L
        for j = 1:N*L
            if abs(G(i,j)) < 0.00001
                G(i,j) = 0;
            end
        end
    end
    
    [U,S,U2] = svds(G,L);
    size(U)
    % VV^T = G. columns are relaxed variables in SDP, V(:,index(i,k))
    V = U(:, 1:L) * sqrt(S(1:L, 1:L));   
    
    basis_vector = zeros(L,1);
    basis_vector(1) = 1;
    rounding = V * (V(0*N+1:0*N+L, 1:L) \ basis_vector)
    
    % print this rounded vector in nice form
    fprintf('\nRounding from the vectors: ');
    for i=1:N
        disp(rounding((i-1)*L+1:(i-1)*L+L)');
    end
    
    Restimated = zeros(N,1);
    % z = zeros(L,1); 
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
        Y(:,i) = circshift(Y(:,i), -Restimated(i));
    end

    z = SDPugJustPhase_pos(Y,sigma,DEBUG);
    
end