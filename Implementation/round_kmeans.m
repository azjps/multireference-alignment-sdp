% V in R^{LN x L} is basis of L-dimensional non-trivial eigenspace of G
% Returns rounded version of V to a single vector which looks close
% to a set of N indicator vectors of length L.
function v_kmeans = round_kmeans(V)
    [N L] = size(V);
    N = N/L;
    I1 = @(i,k) (i-1)*L + k;    % index into 2-D vector v(i,k) in R^{NL}   
    I2 = @(i,k) (k-1)*N + i;    % index into 2-D vector v(k,i) in R^{LN}
    
    best_L1 = 0;
    for rN = 1:N
        rL = 1;
        basis_vector = zeros(L,1);
        basis_vector(rL) = 1;
        v = real(V * (V(I1(rN,1):I1(rN,L), 1:L) \ basis_vector));
        
        % what is L1 norm to set of N indicator vectors of length L?
        
        L1_diff = sum(sign(v));
        for i = 1:N
            M = max(v(I1(i,1):I1(i,L)));
            if M < 1
                L1_diff = L1_diff + 1 - 2*M;
            else
                L1_diff = L1_diff - 1;
            end
        end
        if L1_diff < best_L1 || rN == 1
            v_best = v;
            best_L1 = L1_diff;
        end 
    end
    
    v_kmeans = zeros(N*L, 1);
    X = zeros(N, L);
    for i = 1:N
        X(i, :) = v_best(I1(i,1):I1(i,L))';
    end
    
    [IDX, C] = kmeans(X, L, 'Start', eye(L));    % C is L x N
    for i = 1:N
        v_kmeans(I1(i,IDX(i))) = 1; 
    end
end