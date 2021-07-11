function random_harmonic(L, trials)
    for i = 1:trials
        harmonic = zeros(L,1);
        harmonic(1) = 1;
        harmonic(2:floor(L/2)) = exp(2*pi*1i*rand(L-1,1));
        harmonic(L:floor(L/2)) = 
        disp('New result:');
        disp(harmonic);
        disp(fft(harmonic)/sqrt(L)); 
    end
end