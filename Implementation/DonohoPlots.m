function DonohoPlots

clf

% configurable values
DEBUG = false;           % print status
L = 3;                  % number of rotations
sigmaval = [0:0.1:1];             % noise level
Nval = [5:5:50];      % number of input signals
RUNS = 10;               % number of iterations to run over
ScatS = 100;   % size of scattering - visual

estimator_list = {                              ... 
    @pickfirst,                                 ... % Just pick the first one!
    @angsynchshiftspec,                         ... % Do spectral angular synchornization for the shifts
    @InvarianceMethod,                          ... % Invariance method
    @SDPugJustPhase,                            ... % SDP MLE Synch recovering just shifts
};  


%    @SDPforMLEshiftandsignalMoreConstraints,    ... % SDP MLE Synch recovering signal, using more constraints 
%    @SDPforMLEshiftandsignal,                   ... % SDP MLE Synch recovering signal too
%    @SDPforMLEshiftandsignalNonIndicators,      ... % SDP MLE Synch recovering signal too, non-indicator

METHODS = length(estimator_list);

NvalNum  = length(Nval);
sigmaNum = length(sigmaval);

resultsE = zeros(NvalNum, sigmaNum, METHODS);      % estimator error results
times    = zeros(NvalNum, sigmaNum, METHODS);      % runtime

% returns (positive) complex phase of z, normalized to be between [0,1)
add_one_if_neg = @(z) z + (z < 0);
pos_phase = @(z) add_one_if_neg(atan2(imag(z),real(z))/(2*pi));

    %MakingScattering
         ScatX = zeros(NvalNum*sigmaNum,1);
         ScatY = ScatX;
         ScatC = zeros(NvalNum*sigmaNum,METHODS);
         ScatCounter = 1;
    %MakingScattering


for Ncounter = 1:NvalNum
for sigmacounter = 1:sigmaNum
for run = 1:RUNS
    N = Nval(Ncounter);             
    sigma = sigmaval(sigmacounter);
    
    x = random('normal',0,1,L,1);  % each component with size O(1);
    x_ps     = arrayfun(@(z) abs(z), fft(x)/sqrt(L));
    x_phases = zeros(L,L);
    for l = 0:L-1
        x_phases(:, l+1) = arrayfun(pos_phase, fft(circshift(x,l)));
        % FromFourier(x_ps, x_phases(:, l+1)) % double check
    end
    
    % TODO: print out phases for each rotation of x
    
    % measurements
    Y = zeros(L,N);

    trueShifts = zeros(N,1);
    
    if DEBUG == true
        fprintf('\nRun number: %d, Number of sampled signals: %d\nUnderlying signal: ', ...
            [run; N]);
        disp(x');
        fprintf('Underlying signal Fourier norms, phases: \n');
        disp(x_ps')
        disp(x_phases')
    end
    
    for k = 1:N
        shift = random('unid',L);
        trueShifts(k) = shift;
        Y(:,k) = circshift(x,shift) + sigma*random('normal',0,1,L,1);
        
        if DEBUG == true
             fprintf('Rotation %2d: ', k);
             disp(Y(:,k)');
        end
    end

    for method = 1:METHODS
        % run method and time
        tic
        xest = estimator_list{method}(Y, sigma, DEBUG);
        runtime = toc;
        times(Ncounter,sigmacounter,method) = times(Ncounter,sigmacounter,method) + runtime;
        
        error = distshiftless(x,xest);
        if DEBUG == true
             fprintf('Estimate %d with error %2.3f, time %2.3f: \n', ...
                 [method; error; runtime]);
             disp(xest');
        end
        resultsE(Ncounter,sigmacounter,method) = resultsE(Ncounter,sigmacounter,method) + error;
    end
    
end

clf 

% average and plot
resultsE = resultsE / RUNS;
times = times / RUNS;


   %MakingScattering
         ScatX(ScatCounter) = sigma;
         ScatY(ScatCounter) = N;
         for method = 1:METHODS
            ScatC(ScatCounter,method) = resultsE(Ncounter,sigmacounter,method);
         end
         ScatCounter = ScatCounter +1;

end
end



figure(1)
ScatXB = [];
ScatYB = [];
ScatCB = [];

for method = 1:METHODS
    ScatXB = [ScatXB;ScatX+ceil(sigmaval(sigmaNum))*(method-1)];
    ScatYB = [ScatYB;ScatY];
    ScatCB = [ScatCB;ScatC(:,method)];    
end

scatter(ScatXB,ScatYB,ScatS,ScatCB,'filled')


% for method = 1:METHODS
%     subplot(1,METHODS,method)
%     hold on
%     %MakingScattering
%         scatter(ScatX,ScatY,ScatS,ScatC(:,method),'filled')
%     %MakingScattering
% 
% end


end
   
