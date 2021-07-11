


cvx_begin sdp
    variable G(2, 1)
    variable W(2, 1)
    G(1) == 2
    W(1) == G(1) + G(2)
    W(2) == G(1) - G(2)
    min(sum(G)+sum(W))
cvx_end

G
W