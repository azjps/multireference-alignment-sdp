% I have not tested this because I do not have access to cvx.

% Y is in R^{L x N} is N observations of signal on Z_L 
% N is the number of noisy signals to sample
% sigma is the standard deviation of the additive Gaussian noise
% for example:
%   SDPugJustPhase([1;2;3], 4, 0.5)
function z = SDPugJustPhase(Y, sigma, DEBUG)
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
    A
    
    [U0l,V0,U0r] = svds(A,N+L)
    eig(A) 
    
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
    
    sum(sum(arrayfun(@(z) abs(z), D1)))/L;
    
    % D1(1:N,1:N)
    
    [U1,V1,U7] = svd(D1(1:N,1:N));
    [U2,V2,U8] = svd(D1(N+1:2*N,N+1:2*N));
    [U3,V3,U9] = svd(D1(2*N+1:3*N,2*N+1:3*N));
    
%     A2 = zeros(N*L, N*L);
%     for i = 1:N
%         for j = 1:N
%             for k = 1:L
%                 for l = 1:L
%                     A2(index2(i,k),index2(j,l)) = A(index(i,k),index(j,l));
%                 end
%             end
%         end
%     end
%     inv(dftmtx(N*L)) * A * dftmtx(N*L)
%     inv(dftmtx(N*L)) * A2 * dftmtx(N*L)
      [U,S,U] = svds(A,10);
% %     U
%      S
%     
%     A3 = A;
%     for i = 1:N
%         A3((i-1)*L+1:i*L, (i-1)*L+1:i*L) = eye(L) / L; 
%     end
%     for i = 1:N
%         for j = 1:N 
%             normalization = sum(sum(A3((i-1)*L+1:i*L, (j-1)*L+1:j*L)));
%             for k = 1:L
%                 for l = 1:L
%                     A3(index(i,l),index(j,k)) = A3(index(i,l),index(j,k)) / normalization; % -(norm(y(:,i))^2 + norm(y(:,j))^2)
%                 end
%             end
%         end
%     end
%     A3
%     eig(A3)
%     trace(A3*A)
%     [U,S,U2] = svds(A3,10);
%     S
    

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
%         for i = 1:N
%             for j = i+1:N
%                 for l = 1:L
%                     for k = 1:L
%                         0 <= G(index(i,l),index(j,k)) % <= G(index(i,l),index(i,l))
%                     end
%                 end
%             end
%         end
        % constraint (4)
        % sum G_{ij} = 1
        for i = 1:N
            for j = i+1:N
                sum(sum(G(index(i,1):index(i,L), index(j,1):index(j,L)))) == 1
            end
        end
        
        G >= 0 : Q;
    cvx_end
    
    
    % bunch of verification computations
    G
%     D2 = DFT * G * DFT'
%     [U4,V4,U10] = svd(D2(1:N,1:N));
%     [U5,V5,U11] = svd(D2(N+1:2*N,N+1:2*N));
%     [U6,V6,U12] = svd(D2(2*N+1:3*N,2*N+1:3*N));
%     
%     trace(A*G)
%     trace(D1*D2)
%     trace(D1(1:N,1:N)*D2(1:N,1:N) + D1(N+1:2*N,N+1:2*N)*D2(N+1:2*N,N+1:2*N) + ...
%        D1(2*N+1:3*N,2*N+1:3*N)*D2(2*N+1:3*N,2*N+1:3*N))
%     
%     trace(D1(1:N,1:N)*D2(1:N,1:N))
%     trace(V1(1,1)*V4(1,1)*U1(1:N,1)*U1(1:N,1)'*U4(1:N,1)*U4(1:N,1)')
%     trace(V1(1,1)*V4(1,1)*U4(1:N,1)'*U1(1:N,1)*U1(1:N,1)'*U4(1:N,1))
% 
%     U1(1:N,1)'
%     U4(1:N,1)'
%     
%     U2(1:N,1)'
%     U5(1:N,1)'
%     
%     U3(1:N,1)'
%     U6(1:N,1)'
%     
%     V1(1,1)*V4(1,1)*norm(U4(1:N,1)'*U1(1:N,1))^2 + V2(1,1)*V5(1,1)*norm(U5(1:N,1)'*U2(1:N,1))^2 ... 
%          + V3(1,1)*V6(1,1)*norm(U6(1:N,1)'*U3(1:N,1))^2
%     
%      
%     norm(U4(1:N,1)'*U1(1:N,1))^2
%     norm(U5(1:N,1)'*U2(1:N,1))^2
%     norm(U6(1:N,1)'*U3(1:N,1))^2
%     
%     trace(D1(1:N,1:N))*trace(D2(1:N,1:N))*norm(U4(1:N,1)'*U1(1:N,1))^2 + trace(D1(N+1:2*N,N+1:2*N))*trace(D2(N+1:2*N,N+1:2*N))*norm(U5(1:N,1)'*U2(1:N,1))^2 + ...
%        trace(D1(2*N+1:3*N,2*N+1:3*N))*trace(D2(2*N+1:3*N,2*N+1:3*N))*norm(U6(1:N,1)'*U3(1:N,1))^2
      
   
   
%     [U,S,U2] = svds(D(1:N,1:N));
%     U
%     S
%     [U,S,U2] = svds(D(N+1:2*N,N+1:2*N));
%     U
%     S
%     [U,S,U2] = svds(D(2*N+1:3*N,2*N+1:3*N));
%     U
%     S
    
    
%     G2 = zeros(N*L, N*L);
%     for i = 1:N
%         for j = 1:N
%             for k = 1:L
%                 for l = 1:L
%                     G2(index2(i,k),index2(j,l)) = G(index(i,k),index(j,l));
%                 end
%             end
%         end
%     end
%     inv(dftmtx(N*L)) * G * dftmtx(N*L)
%     inv(dftmtx(N*L)) * G2 * dftmtx(N*L)
     [U,S,U2] = svd(G);
%     U
%     S
    
    for i = 1:N*L
        for j = 1:N*L
            if abs(G(i,j)) < 0.00001
                G(i,j) = 0;
            end
        end
    end
    
    [U,S,U2] = svds(G);
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