% Y in R^{L x N}
function z = KmeansCluster(Y, sigma, DEBUG)
    [L N] = size(Y);
    I1 = @(i,k) (i-1)*L + k;    % index into 2-D vector v(i,k) in R^{NL}   
    % I2 = @(i,k) (k-1)*N + i;    % index into 2-D vector v(k,i) in R^{LN}
    
    %% populate data
    X = zeros(N*L, L);
    for i = 1:N
        for k = 1:L
            X(I1(i,k), :) = circshift(Y(:,i), k); 
        end
    end
    
    [IDX, C] = kmeans(X, L, 'Start', X(1:L, :));    % C is L x N
    C
    z = C(:, 1);
    % zeros(L, 1);
    % count = 0;
    % for i = 1:N*L
    %     if IDX(i) == 1
    %         z = z + X(i)';
    %         count = count + 1;
    %     end
    % end
    z = real(PowerSpectrumCorrection(z, Y, sigma));
end