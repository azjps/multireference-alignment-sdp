% use both objective function with just rotations and 
% objective function with Fourier phase of signal
% combined

function z = SDPforMLEshiftandsignalNonIndicators(Y,sigma,DEBUG)

lambda = 0; % mixing constant

[L N] = size(Y);

z = random('normal',0,1,L,1);  
z = PowerSpectrumCorrection(z,Y,sigma); 

Fz = fft(z);
Amod = abs(Fz) / sqrt(L);       % power spectrum

%COMPUTING LIKELIHOODS

% SOLVING THE SDP
% index into the big matrix
% v and w are 0/1 indicator variables.
% v is 1 if index refers to an indicator of shift: 
%   (i,k) = i-th signal has k-th rotation  
% w is 1 if index refers to a Fourier phase
%   2*j-1 = real part of jth Fourier phase, 
%   2*j   = imag part of jth Fourier phase.
index = @(v,w,i,k,j) v*((i-1)*L + k) + w*(L*N+j); 
                        
C = zeros(N*L+2*L,N*L+2*L);

for i = 1:N         % signal #
    Fy = fft(Y(:,i))/sqrt(L);
    for l = 1:L     % rotation 
        ind1 = index(1,0,i,l,0);
        for h = 1:L % hth component of Fourier phase of x
            % Compute likelihood of R(-l)y_i and the h'th harmonic of x

            ab = Fy(h)*exp(2*pi*1i*(l-1)*(h-1)/L);
            a = real(ab);
            b = imag(ab);

            ind2_re = index(0,1,0,0,2*h-1);
            ind2_im = index(0,1,0,0,2*h);

            C(ind1,ind2_re) = a;
            C(ind2_re,ind1) = a;
            C(ind1,ind2_im) = b;
            C(ind2_im,ind1) = b;
        end
    end
end

if lambda ~= 0
    D = zeros(N,N,L);
    for i = 1:N
        for j = 1:N     % or use symmetry with j = i+1:N
            for l = 1:L
                D(i,j,l) = circshift(Y(:,i), -l+1)' * Y(:,j);
            end
        end
    end

    % We can relax to a semidefinite program as follows. 
    % Minimizing sum_ij D(i,j) = sum_ij sum_k Dhat(i,j)G(@ik,@jk)/L
    % is similar to the SDP:
    % min_{l_1, ..., l_n in [L]} tr(AG)

    for i = 1:N
        for j = i+1:N
            for k = 1:L
                for l = 1:L
                    C(index(1,0,i,l,0),index(1,0,j,k,0)) = 2*lambda*D(i,j,mod(l-k,L)+1); % -(norm(y(:,i))^2 + norm(y(:,j))^2)
                    C(index(1,0,j,k,0),index(1,0,i,l,0)) = C(index(1,0,i,l,0),index(1,0,j,k,0));
                end
            end
        end
    end
end

C

cvx_begin sdp
    variable G(N*L+2*L,N*L+2*L) symmetric
    dual variable Q
    %We want to maximize Tr(CG) over PSD symmetric G
    maximize( trace(C*G) )
    
    for i=1:N
        % A is the L*L block identity matrix in the diagonal of a 
        % N*L + 2*L square matrix, used to constrain diagonal elements
        A = zeros(N*L+2*L,N*L+2*L);
        for k1=1:L
            % first N*L diagonal elements equal to 1/L
            ind1 = index(1,0,i,k1,0);
            % remove ambiguity from rotations
            % if i == 1 && k1 == 1
            %     G(ind1,ind1) == 1
            % elseif i == 1
            %     G(ind1,ind1) == 0
            % else
            %     A(ind1,ind1) = 1;
                G(ind1,ind1) == 1/L
            % end

            % Within each i, non diagonal entries need to be zero.
            for k2=(k1+1):L
                ind2 = index(1,0,i,k2,0);
                G(ind1,ind2) == 0
                G(ind2,ind1) == 0
            end
        end
        
        % if i ~= 1
        %     trace(A*G) == 1
        % end
    end

    for j=1:L
        %The Power Spectrum of Z
        inda = index(0,1,0,0,2*j-1);
        indb = index(0,1,0,0,2*j);
        G(inda,inda)+G(indb,indb) == Amod(j)
    end
    
    G >= 0: Q; 
cvx_end

G

% rank-1 approximation?
% [U,S,U2] = svds(G);
% U
% S
% G(L*N+1:(N+2)*L, L*N+1:(N+2)*L)
% [U,S,U2] = svds(G(L*N+1:(N+2)*L, L*N+1:(N+2)*L));
% U
% S

% playing around with XX^*
% cov_X1 = G(L*N+1:(N+2)*L, L*N+1:(N+2)*L)
% cov_X2 = zeros(L,L);
% for i = 1:L
%     for j = 1:L
%         cov_X2(i,j) = (cov_X1(2*i-1,2*j-1)+cov_X1(2*i,2*j)) ...
%             + (cov_X1(2*i,2*j-1)-cov_X1(2*i-1,2*j)) * 1i;
%     end
% end
% cov_X2


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% random rounding approach
if DEBUG == true
    [U,S,U2] = svd(G);      % G = U*S*V'
    V = U*sqrt(S);          % rows of V are the vectors!
    g = random('normal',0,1,N*L+L*2,1);
    rounding = V * g;

    Restimated = zeros(N,1);
    z = zeros(L,1);
    for i = 1:N
        best = -1000;
        for k = 1:L
            prod = rounding(index(1,0,i,k,0));
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
    fprintf('\nRecovered signal from averaging of rotations: \n')
    disp(z')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for i = 1:L
G2 = G(L*N+1:L*N+2*L,L*N+1:L*N+2*L);
[U,S,U2] = svds(G2);
V = U*sqrt(S);   % VV^T = G. columns are relaxed variables in SDP, V(:,index(i,k))
sizeV = size(V);

add_one_if_neg = @(z) z + (z < 0);
pos_phase = @(z) add_one_if_neg(atan2(imag(z),real(z))/(2*pi));
phi = zeros(L,1);

% random rounding approach
g = normrnd(0,1,sizeV(2),1);
rounding = V * g;

% Note: there is a sign ambiguity when recovering the phases
% We use the DC component of the original Y
% and compare it against the DC component of the recovered Fourier
% coefficients X to determine the sign
DC_Y = sum(sum(Y));
DC_sign = 1;
if (DC_Y > 0 && rounding(1) < 0) || (DC_Y < 0 && rounding(1) > 0)
    DC_sign = -1;
end    
    
z = zeros(L,1);
for j = 1:L
    z(j) = DC_sign * (rounding(2*j-1) + 1i*rounding(2*j));
    phi(j) = pos_phase(z(j));
end

% z = PowerSpectrumCorrection(z,Y,sigma); 
z = FromFourier(Amod, phi);
if DEBUG == true
    fprintf('Recovered signal from phases, and phases: \n')
    disp(z')
    disp(phi')
end   

% end
end
  
