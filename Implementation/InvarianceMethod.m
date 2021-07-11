% Y is in R^{L x N} is N observations of signal on Z_L 
% N is the number of noisy signals to sample
% sigma is the standard deviation of the additive Gaussian noise
function min_z = InvarianceMethod(Y,sigma, DEBUG)

[L N] = size(Y);  

D = 1000;           % discretization of first Fourier phase
MAX = ceil(D/L);    % number of discretized units to run to

z = random('normal',0,1,L,1);               % REMOVE REMOVE REMOVE
z = PowerSpectrumCorrection(z,Y,sigma);
Fz = fft(z);
Amod = abs(Fz) / sqrt(L);         % power spectrum

FY = fft(Y);

Z = zeros(MAX,L);

DCcomp = FY(1,:) * ones(N,1) / N;   % average of DC components

min_dist = -1;
min_z = zeros(L,1);
for j = 1:MAX
    Z(j,0+1) = Amod(0+1) * (DCcomp/abs(DCcomp));
    Z(j,1+1) = Amod(1+1) * exp(2*pi*1i*j/500);

    Inv = zeros(L,1);

    for k=2:L-1
        InvY = (FY(1+1,:).^(-k)).*FY(k+1,:);
        Inv(k) = InvY*ones(N,1)/N;
        Z(j,k+1) = Z(j,1+1)^k*Inv(k);
        Z(j,k+1) = (Z(j,k+1)/abs(Z(j,k+1)))*Amod(k+1);    
    end
    
    z = ifft(Z(j, :))'*sqrt(L);
    z = PowerSpectrumCorrection(z,Y,sigma);
    new_dist = distshiftless(z,Y);
    % fprintf('Distance at discretization %f: %f\n', [j/D; new_dist]);
    
    if j == 1 || new_dist < min_dist
        min_dist = new_dist;
        min_z = z;
    end
end

end
  
