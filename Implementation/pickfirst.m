function z = pickfirst(Y,sigma,DEBUG)
    z = PowerSpectrumCorrection(Y(:,1), Y, sigma);
end