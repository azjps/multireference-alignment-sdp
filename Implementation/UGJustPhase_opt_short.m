% Y is in R^{L x N} is N observations of signal on Z_L 
% sigma is the standard deviation of the additive Gaussian noise
% DEBUG is a boolean flag as to whether to print debug information
function z = UGJustPhase_opt_short(Y, sigma, DEBUG)
    [L N] = size(Y);
    I1 = @(i,k) (i-1)*L + k;    % index into 2-D vector v(i,k) in R^{NL}   
    I2 = @(i,k) (k-1)*N + i;    % index into 2-D vector v(k,i) in R^{LN}
    P = zeros(N*L,N*L);         % permutation matrix
    for i = 1:N
        for k = 1:L
            P(I1(i,k), I2(i,k)) = 1;
        end
    end
    DFT = kron(conj(dftmtx(L)) / sqrt(L), eye(N)) * P';     % block DFT

    %% componentwise fft 
    V_DFT = zeros(N*L,L);       % C = V_DFT*V_DFT'/L is data Gram matrix
    for i = 1:N                 % gram(k) = <R{k-1}*Y_i, Y_1>
        gram = arrayfun(@(k) circshift(Y(:,i), k-1)' * Y(:,1), 1:L); 
        Fc = fft(gram);
        for k = 1:L
            V_DFT(I2(i,k), k) = sign(Fc(k));
        end
    end
    
    %% rounding one signal's vector to indicator vector
    V = DFT' * V_DFT
    [rL, rN] = deal(1,1);
    basis_vector = zeros(L,1);
    basis_vector(rL) = 1;
    v = real(V * (V(I1(rN,1):I1(rN,L), 1:L) \ basis_vector));
    
    %% choose best rotation for each image
    z = zeros(L,1); 
    for i = 1:N
        [M,R] = max(v(I1(i,1):I1(i,L)));
        z = z + circshift(Y(:,i), -R+1);
    end
    z = real(PowerSpectrumCorrection(z/N,Y,sigma));
end