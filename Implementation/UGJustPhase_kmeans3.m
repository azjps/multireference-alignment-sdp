function z = UGJustPhase_kmeans3(Y, sigma, DEBUG)
    [L N] = size(Y);            % Y in R^{L x N} is N observations of Z_L signal  
    I = @(i,k) (i-1)*L + k;     % index into 2-D vector v(i,k) in R^{NL}   
   
    %% componentwise fft 
    V = zeros(N*L,L);
    for i = 1:N                 % gram(k) = <R{k-1}*Y_i, Y_1>
        gram = arrayfun(@(k) circshift(Y(:,i), -k)' * Y(:, 1), 0:L-1); 
        V(I(i,1):I(i,L),1) = fft(sign(ifft(gram)))' / L;
        for k = 2:L
            V(I(i,1):I(i,L),k) = circshift(V(I(i,1):I(i,L),1), k-1);
        end
    end
    
    %% rounding: pick vector from eigenspace basis V resembling indicator
    v = round_kmeans(V);        % naive rounding: v = real(V*ones(L,1));
    z = zeros(L,1); 
    for i = 1:N
        [M,R] = max( v(I(i,1):I(i,L)) );
        z = z + circshift(Y(:,i), -R+1);
    end
    z = real(PowerSpectrumCorrection(z/N, Y, sigma));
end