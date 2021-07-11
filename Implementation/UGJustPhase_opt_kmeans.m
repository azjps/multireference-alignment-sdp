% Y in R^{L x N} is N observations of signal on Z_L 
% sigma is the standard deviation of the additive Gaussian noise
% DEBUG is a boolean flag as to whether to print debug information
function z = UGJustPhase_opt_kmeans(Y, sigma, DEBUG)
    TIMING = false;
    [L N] = size(Y);
    I1 = @(i,k) (i-1)*L + k;    % index into 2-D vector v(i,k) in R^{NL}   
    I2 = @(i,k) (k-1)*N + i;    % index into 2-D vector v(k,i) in R^{LN}
    P = zeros(N*L,N*L);         % permutation matrix
    for i = 1:N
        for k = 1:L
            P(I1(i,k), I2(i,k)) = 1;
        end
    end
    DFT = kron(dftmtx(L)' / sqrt(L), eye(N)) * P'     % block DFT
    kron(eye(N), dftmtx(L)' / sqrt(L))
    
    %% componentwise fft 
    if TIMING
        tic
    end
    g_DFT = zeros(N*L,L);
    V_DFT = zeros(N*L,L);       % C = V_DFT*V_DFT'/L is data Gram matrix
    for i = 1:N                 % gram(k) = <R{k-1}*Y_i, Y_1>
        gram = arrayfun(@(k) circshift(Y(:,i), k-1)' * Y(:,1), 1:L); 
        g_DFT(I1(i,1):I1(i,L), :) = ones(L,1) * sign(fft(gram));
        c = fft(gram);
        for k = 1:L
            V_DFT(I2(i,k), k) = sign(c(k));
        end
    end
    if TIMING
        fprintf('Time for fft: %f\n', toc);
    end
    
    %% rounding one signal's vector to indicator vector
    if TIMING
        tic
    end
    
    V_DFT
    V = DFT' * V_DFT
    V - kron(ones(N,1), dftmtx(L)/sqrt(L)) .* g_DFT
    V * ones(L,1)
    % round_one, round_random, round_best, round_bestL1
    v = round_kmeans(V);
    if TIMING
        fprintf('Time for rounding: %f\n', toc);
    end
    
    %% choose best rotation for each image
    if TIMING
        tic
    end
    z = zeros(L,1); 
    for i = 1:N
        [M,R] = max(v(I1(i,1):I1(i,L)));
        z = z + circshift(Y(:,i), -R+1);
    end
    z = real(PowerSpectrumCorrection(z/N, Y, sigma));
    if TIMING
        fprintf('Time for argmax: %f\n', toc);
    end
end