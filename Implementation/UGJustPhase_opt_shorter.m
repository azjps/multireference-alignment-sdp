function z = UGJustPhase_opt_shorter(Y, sigma, DEBUG)
    [L N] = size(Y);            % Y in R^{L x N} is N observations of Z_L signal  
    I = @(i,k) (i-1)*L + k;     % index into 2-D vector v(i,k) in R^{NL}   
   
    %% componentwise fft 
    g_DFT = zeros(N*L,L);
    for i = 1:N                 % gram(k) = <R{k-1}*Y_i, Y_1>
        gram = arrayfun(@(k) circshift(Y(:,i),k)' * Y(:,1), 0:L-1); 
        g_DFT(I(i,1):I(i,L), :) = ones(L,1) * sign(fft(gram));
    end
    V = kron(ones(N,1), dftmtx(L) / sqrt(L)) .* g_DFT;
    
    %% rounding: pick vector from eigenspace basis V resembling indicator
    v = round_bestL1(V);        % naive rounding: v = real(V*ones(L,1));
    z = zeros(L,1); 
    for i = 1:N
        [M,R] = max( v(I(i,1):I(i,L)) );
        z = z + circshift(Y(:,i), -R+1);
    end
    z = real(PowerSpectrumCorrection(z/N, Y, sigma));
end