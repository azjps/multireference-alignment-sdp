% Y in R^{L x N}
% returns what the actual shifts / fourier phases should be,
% so that it may be tested with the objective function
function [u,phi] = ActualShifts(Y)
	x = Y(:,1);
    u = zeros(N,1);
	
	for n = 1:N
        yn = Y(:,n);
        min_dist = norm(x - yn);
	    for l = 2:L
            new_dist = norm(x - circshift(yn, l-1));
            if min_dist > new_dist
                min_dist = new_dist;
                u(n) = l-1;
            end
		end
	end
	
    add_one_if_neg = @(z) z + (z < 0);
    pos_phase = @(z) add_one_if_neg(atan2(imag(z),real(z))/(2*pi));
    
    phi = arrayfun(pos_phase, fft(x)/sqrt(L));
end