
% Compute the first 2*T risk neutral probabilities
% for which the resulting yield curve is flat at r0.
% Initial interest rate is r0, and each time step
% of 1/2 years, either goes up by factor of u
% or goes down by factor of d.

function Q = SemiAnnualFlatQ(r0, u, d, T)
    asset1 = @(r) 1/(1+(u*r)/2) - 1/(1+(d*r)/2);
    asset2 = @(r) 1/(1+(d*r)/2);

    R = zeros(2*T+1, 2*T+1);
    R(1,1) = r0;        % upper triangular binomial interest rate tree
    for t = 2 : 2*T+1
        R(1,t) = R(1,t-1) * u;
    end
    for t = 2 : 2*T+1
        for s = 1 : t-1
            R(s+1,t) = R(s,t) * (d/u);
        end
    end
    
    Q = zeros(2*T-1, 1);
    for t = 1/2 : 1/2 : T
        t2 = 2*t;
        
        P  = (1 + r0/2)^(-2*t);
        P1 = SemiAnnualInterestDeriv(asset1, R, Q(1:t2-1,1));
        P2 = SemiAnnualInterestDeriv(asset1, R, Q(1:t2-1,1));
        Q(t2) = (P - P2) / P1;
    end
end