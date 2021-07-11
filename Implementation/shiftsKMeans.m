% V = [v_1 v_2 ... v_L]

function [] = shiftsKMeans(V,eps)

L = length(V);

% minimize L1 norm
cvx_begin
    variable x(L);
    minimize( norm(V*x,1) );
    subject to
        V(1,:)*x == 1;
cvx_end

Vopt = V*x;

% k-means cluster

% cluster mean

end