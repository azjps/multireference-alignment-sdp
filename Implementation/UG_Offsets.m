function Z = UG_Offsets(X,Y,sigma,DEBUG)
	[L N] = size(Y);
    D = zeros(N,N,L);
    for i = 1:N
        for j = 1:N     % or use symmetry with j = i+1:N
            for l = 1:L
                D(i,j,l) = norm(circshift(Y(:,i), -l+1) - Y(:,j))^2;
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
        for j = i+1:N
            for k = 1:L
                for l = 1:L
                    A(index(i,l),index(j,k)) = D(i,j,mod(l-k,L)+1); % -(norm(y(:,i))^2 + norm(y(:,j))^2)
                end
            end
        end
    end
    
    % D
    % A
    
    cvx_begin sdp
        variable G(N*L, N*L) symmetric
        dual variable Q
        minimize( trace(G*A) )
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
    
    Z = zeros(N,N);
    for n1 = 1:N
        for n2 = 1:N
            % entry l corresponds to rotation (l-1)
            likelihood = arrayfun(@(l) G(index(n1,1+l), index(n2,1)), 0:L-1);
            [maxlikelihood, k] = max(likelihood);
            Z(n1,n2) = mod(k-1,L);
        end
    end
end