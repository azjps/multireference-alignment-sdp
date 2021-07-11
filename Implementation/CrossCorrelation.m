
function z = CrossCorrelation(Y, sigma, DEBUG)
    [L N] = size(Y); 
    %% shifts of y_1
    Y1 = zeros(L, L);
    for l = 1:L
        Y1(:,l) = circshift(Y(:,1), l-1);
    end
    
    %% 
    z = Y1(:,1);
    Restimated = zeros(N, 1);
    for i = 2:N
        best = -1000;
        for k = 1:L
            prod = Y1(:,k)' * Y(:,i);
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