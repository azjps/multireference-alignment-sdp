function z = angsynchshiftspec(Y,sigma,DEBUG)

[L N] = size(Y);

%z = ones(L,1);

RelShifts = 73*ones(N,N); 

for i=1:(N-1)
for j=(i+1):N    
     %Y(:,i)
     %Y(:,j)
     optshift = optshiftfrom2to1(Y(:,i),Y(:,j));
     RelShifts(i,j) = optshift;
     RelShifts(j,i) = -RelShifts(i,j);
end
end

W = exp(2*pi*1i*RelShifts/L);
W = W + N*eye(N);

[g lambda] = eigs(W,1);

for i=1:N
   g(i) = g(i)/norm(g(i));
   g(i) = angle(g(i));
   g(i) = L*g(i)/(2*pi);
   g(i) = round(g(i));
   g(i) = mod(g(i),L);
   yc = Y(:,i);
   Y(:,i) = circshift(yc,-g(i));
end



z = Y*ones(N,1)/N;

z = PowerSpectrumCorrection(z,Y,sigma);



end