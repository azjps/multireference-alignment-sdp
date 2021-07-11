function z = SDPforMLEshiftandsignal(Y,sigma, DEBUG)

[L N] = size(Y);

z = random('normal',0,1,L,1);  %REMOVE REMOVE REMOVE
z = PowerSpectrumCorrection(z,Y,sigma); %NOW IT HAS THE CORRECT NORM

Fz = 1/sqrt(L)*fft(z);
Amod = abs(Fz);

M = L;

%COMPUTING LIKELIHOODS





%SOLVING THE SDP
    index = @(v,w,i,k,j,m) v*((i-1)*L + k) + w*(L*N+(j-1)*M+m); 
                % index into the big matrix
                % v and w are indicators and the rest is the index
                % First there are the variables for the shifts and then the
                % ones for the phases of the fourier transform of x

                
                
                
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
           Xharmonic = 1/sqrt(L)*Amod(j)*exp(2*pi*1i*m/M)*exp(2*pi*1i*(j-1)*[0:L-1]'/L);
           ind2 = index(0,1,i,k,j,m);
           C(ind1,ind2) = real(beta'*Xharmonic);
           indicesRead(ind1,ind2) = 1;
           C(ind2,ind1) = C(ind1,ind2);
           indicesRead(ind2,ind1) = indicesRead(ind1,ind2);
        end
      end
    end
  end

  

%Constraints:

%rest M*L diagonal elements equal to 1/M


 
    cvx_begin sdp
        variable G(N*L+L*M,N*L+L*M) symmetric
        dual variable Q
        %We want to maximize Tr(CG) over PSD simetric G
        maximize( trace(C*G) )

        for i=1:N
        for k1=1:L
            %first N*L diagonal elements equal to 1/L
            ind1 = index(1,0,i,k1,0,0);
            G(ind1,ind1) == 1/L
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
% 		for i = 1:N*L+L*M,N*L+L*M
%             for j = i+1:N*L+L*M,N*L+L*M
%                 0 <= G(i,j) % <= G(index(i,l),index(i,l))
%             end
%         end
        
        
         G >= 0: Q; 
        %This is just PSD!
       
    cvx_end

    indicesRead;
    
%G

[U,S,V] = svd(G);   % G = U*S*V'
% rows of V are the vectors!!!




%ROUNDING THE SOLUTION  - NAIVELY

 %pick random g gaussian with size N*L+L*M
    g = random('normal',0,1,N*L+L*M,1);
 %only need to round the signal part
    %for each j pick the m that agrees more with g
    
    %The vectors are the rows of V
    
  phi = zeros(L,1);  %will contain the phases 0,1/M,...,1
 
    for j=1:L
            m=M;
            phi(j) = 1;
            oldprob = V(index(0,1,0,0,j,m),:)*g; %testing first the phase 0
       for m=1:M-1
            newprob = V(index(0,1,0,0,j,m),:)*g;  %testing every phase
            if newprob>oldprob
                phi(j) = m/M;
                oldprog = newprob;
            end
       end
    end

    %Building z by power spectrum recovery
    z = zeros(L,1);
    Fz = zeros(L,1);
    for j = 1:L
        Fz(j) = Amod(j)*exp(2*pi*1i*phi(j));
        z = z + 1/sqrt(L)*Fz(j)*exp(2*pi*1i*(j-1)*[0:L-1]'/L);
    end
    
%    MustBeZero = norm(Fz - 1/sqrt(L)*fft(z))
        
       %Corz = PowerSpectrumCorrection(z,Y,sigma); %GETTING THE MODULUS RIGHT!

           %MustBeZero = norm(z - Corz);

       
 %     QuotientsOfTheNorms = abs(fft(z)./fft(Corz))
       
    
end
  
