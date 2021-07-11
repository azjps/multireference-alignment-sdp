% V in R^{LN x L} is basis of L-dimensional non-trivial eigenspace of G
% Returns rounded version of V to a single vector which looks close
% to a set of N indicator vectors of length L.
function v = round_random(V)
    [N L] = size(V);
    N = N/L;
    I1 = @(i,k) (i-1)*L + k;    % index into 2-D vector v(i,k) in R^{NL}   
    I2 = @(i,k) (k-1)*N + i;    % index into 2-D vector v(k,i) in R^{LN}
    
    [rL, rN] = deal(1, ceil(N*rand(1)));
    basis_vector = zeros(L,1);
    basis_vector(rL) = 1;
    v = real(V * (V(I1(rN,1):I1(rN,L), 1:L) \ basis_vector));
end