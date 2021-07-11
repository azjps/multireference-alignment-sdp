% recover a signal from its sequence of Fourier amplitudes (power spectrum)
% and phases

% NOTE: THESE ARE INDEXED STARTING AT 0!
% Hence amps(1) corresponds to the amplitude of the 0th Fourier
% coefficient, etc
function x = FromFourier(amps, phases) 
    L = length(amps);

    x = zeros(L,1);
    for n = 1:L
        for k = 1:L
            x(n) = x(n) + amps(k) * exp(2*pi*1i*(phases(k)+(k-1)*(n-1)/L));
        end
        x(n) = x(n) / sqrt(L);  % normalization
    end

end