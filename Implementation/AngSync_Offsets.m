function Z = AngSync_Offsets(X,Y,sigma,DEBUG)
    [L, N] = size(Y);
    Z = zeros(N,N);
    for n1 = 1:N
        for n2 = 1:N
            % entry l corresponds to rotation (l-1)
            likelihood = arrayfun(@(l) norm(Y(:,n1)-circshift(Y(:,n2),l-1)), 1:L);
            [maxlikelihood, k] = min(likelihood);
            Z(n1,n2) = k-1;
        end
    end
end