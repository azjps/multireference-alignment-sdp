% Y is in R^{L x N} is N observations of signal on Z_L 
% N is the number of noisy signals to sample
% sigma is the standard deviation of the additive Gaussian noise
function z = SDPforMLEshiftandsignal(Y,sigma, DEBUG)

[L N] = size(Y);  

z = random('normal',0,1,L,1);               % REMOVE REMOVE REMOVE
z = PowerSpectrumCorrection(z,Y,sigma);     % NOW IT HAS THE CORRECT NORM

Fz = fft(z);
Amod = abs(Fz);         % power spectrum

M = 10;      % discretization of Fourier phases of signal

% SOLVING THE SDP
% index into the big matrix
% v and w are 0/1 indicator variables.
% v is 1 if index refers to an indicator of shift: 
%   (i,k) = i-th signal has k-th rotation  
% w is 1 if index refers to an indicator of a discretization of a Fourier
% phase:
%   (j,m) = j-th Fourier phase has mth discrete value on unit circle 
index = @(v,w,i,k,j,m) v*((i-1)*L + k) + w*(L*N+(j-1)*M+m); 
                
C = zeros(N*L+L*M,N*L+L*M);

indicesRead = C;

for i=1:N
    for k=1:L
        beta = Y(:,i);
        beta = circshift(beta,-k);
        ind1 = index(1,0,i,k,0,0);
        for j=1:L
            for m=1:M
                % Compute likelihood of R(-k)y_i and the j'th harmonic
                % of x with value phi(j) = m/M
                Xharmonic = Amod(j)*exp(2*pi*1i*(m-1)/M)*exp(2*pi*1i*(j-1)*[0:L-1]'/L);
                ind2 = index(0,1,i,k,j,m);
                C(ind1,ind2) = real(beta'*Xharmonic);
                indicesRead(ind1,ind2) = 1;
                C(ind2,ind1) = C(ind1,ind2);
                indicesRead(ind2,ind1) = indicesRead(ind1,ind2);
            end
        end
    end
end

% C  

%Constraints:

% rest M*L diagonal elements equal to 1/M

cvx_begin sdp
    variable G(N*L+L*M,N*L+L*M) symmetric
    dual variable Q
    %We want to maximize Tr(CG) over PSD simetric G
    maximize( trace(C*G) )

    for i=1:N
        A = zeros(N*L+M*L,N*L+M*L);
        for k1=1:L
            %first N*L diagonal elements equal to 1/L
            ind1 = index(1,0,i,k1,0,0);
            
            % remove ambiguity from rotations
            if i == 1 && k1 == 1
                G(ind1,ind1) == 1
            elseif i == 1
                G(ind1,ind1) == 0
            else
                A(ind1,ind1) = 1;
            %    G(ind1,ind1) == 1/L
            end
            
            indicesRead(ind1,ind1) = 30;
        %Within each i, non diagonal entries need to be zero.
        for k2=(k1+1):L
            ind2 = index(1,0,i,k2,0,0);
            G(ind1,ind2) == 0
            G(ind2,ind1) == 0
            indicesRead(ind1,ind2) = 20;
            indicesRead(ind2,ind1) = 20;
        end
        end
        
        if i ~= 1
            trace(A*G) == 1
        end
    end

    for j=1:L
        for m1=1:M
            %second L*M diagonal elements equal to 1/M
            ind1 = index(0,1,0,0,j,m1);
            G(ind1,ind1) == 1/M
            indicesRead(ind1,ind1) = 30;
        %Within each  j, non diagonal entries need to be zero.
        for m2=(m1+1):M
            ind2 = index(0,1,0,0,j,m2); 
            G(ind1,ind2) == 0
            G(ind2,ind1) == 0
            indicesRead(ind1,ind2) = 20;
            indicesRead(ind2,ind1) = 20;
        end
        end
    end

    %Making the matrix non-negative

        % ANDY, CAN YOU ADD THE CONSTRAINT THAT THE ENTRIES OF G ARE
        % NON-NEGATIVE? I THINK THAT IS THE ONLY THING MISSING AND I AM
        % NOT ABLE TO DO IT FOR SOME PROBLEM... JUST ADD THAT AND RUN
        % THE ShiftLessRecovery.m       THANKS!
        % PS: THE PDF ALSO HAS EXPLANATION OF THE METHOD 

    % Non-negativity constraint the lazy way,
    % have not tested
%    for i = 1:N*L+L*M,N*L+L*M
%        for j = i+1:N*L+L*M,N*L+L*M
%            0 <= G(i,j) % <= G(index(i,l),index(i,l))
%        end
%    end


     G >= 0: Q; 
cvx_end

G;

[U,S,U2] = svd(G);      % G = U*S*V'
V = U*sqrt(S);          % rows of V are the vectors!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% random rounding approach
g = random('normal',0,1,N*L+L*M,1);
rounding = V * g
Restimated = zeros(N,1);
z = zeros(L,1);
for i = 1:N
    best = -1000;
    for k = 1:L
        prod = rounding(index(1,0,i,k,0,0));
        if prod > best
            Restimated(i) = k - 1;
            best = prod;
        end
    end
    z = z + circshift(Y(:,i), -Restimated(i));
end
z = z / N;

% random rank-1 decomposition
% g = normrnd(0,1,N*L,1);
% R'
% (G*g)'

z = PowerSpectrumCorrection(z,Y,sigma);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ROUNDING THE SOLUTION  - NAIVELY

% pick random g gaussian with size N*L+L*M
g = random('normal',0,1,N*L+L*M,1);
% Only need to round the signal part
% for each j pick the m that agrees more with g
% The vectors are the rows of V
    
phi = zeros(L,1);  % will contain the phases 0, 1/M, ..., 1

for j = 1:L
    m = M;
    phi(j) = 1;
    oldprob = rounding(index(0,1,0,0,j,m)); % testing first the phase 0
    for m = 1:M-1
        newprob = rounding(index(0,1,0,0,j,m)); % testing every phase
        if newprob > oldprob
            phi(j) = m/M;
            oldprob = newprob;
        end
    end
end

phi

%Building z by power spectrum recovery
z = zeros(L,1);
Fz = zeros(L,1);
for j = 1:L
    Fz(j) = Amod(j)*exp(2*pi*1i*phi(j));
    z = z + 1/sqrt(L)*Fz(j)*exp(2*pi*1i*(j-1)*[0:L-1]'/L);
end
    
end
  
