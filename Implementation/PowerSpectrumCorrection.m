function z = PowerSpectrumCorrection(z,Y,sigma)

    % THIS IS NOT WORKING!!!

    [L N] = size(Y);

    FY = fft(Y);
    FY = abs(FY);

    modF = FY*ones(N,1)/N;

    modF = sqrt(modF.^2 - sigma^2);

    Fz = fft(z);

    for k=1:L
        Fz(k) = (Fz(k)/abs(Fz(k))) * modF(k);
    end

    z = ifft(Fz);

end