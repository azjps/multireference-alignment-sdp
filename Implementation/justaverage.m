function z = justaverage(Y,sigma)

[L N] = size(Y);

z = Y*ones(N,1)/N;

z = PowerSpectrumCorrection(z,Y,sigma);

end