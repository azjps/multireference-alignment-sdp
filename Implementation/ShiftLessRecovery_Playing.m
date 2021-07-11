function ShiftLessRecovery

clf

% configurable values
DEBUG = true;           % print status
L = 2;                  % number of rotations
sigma = 0.5;             % noise level
Nval = 1*[10:10:50];      % number of input signals
RUNS = 50;               % number of iterations to run over

estimator_list = {                          ... 
    @pickfirst,                             ... % Just pick the first one!
    @angsynchshiftspec,                     ... % Do spectral angular synchornization for the shifts
    @InvarianceMethod,                      ... % Invariance method
    @SDPugJustPhase,                        ... % SDP MLE Synch recovering just shifts
    @SDPforMLEshiftandsignalNonIndicators,  ... % SDP MLE Synch recovering signal too, non-indicator
};
%    @SDPforMLEshiftandsignal,               ... % SDP MLE Synch recovering signal too
METHODS = length(estimator_list);

NvalNum  = length(Nval);
resultsE = zeros(NvalNum, METHODS);      % estimator error results
times    = zeros(NvalNum, METHODS);      % runtime

% returns (positive) complex phase of z, normalized to be between [0,1)
add_one_if_neg = @(z) z + (z < 0);
pos_phase = @(z) add_one_if_neg(atan2(imag(z),real(z))/(2*pi));

for Ncounter = 1:NvalNum
for run = 1:RUNS
    N = Nval(Ncounter);             
    
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
        times(Ncounter,method) = times(Ncounter,method) + runtime;
        
        error = distshiftless(x,xest);
        if DEBUG == true
             fprintf('Estimate %d with error %2.3f, time %2.3f: \n', ...
                 [method; error; runtime]);
             disp(xest');
        end
        resultsE(Ncounter,method) = resultsE(Ncounter,method) + error;
    end
end
end
   
% average and plot
resultsE = resultsE / RUNS;
times = times / RUNS;

figure(73);
COLORS = 'rygcbmk';

subplot(1,2,1)
hold on

for method = 1:METHODS
    plot(Nval, resultsE(:,method), COLORS(method))
end
title('Error varying with N')

subplot(1,2,2)
hold on
for method = 1:METHODS
    plot(Nval, times(:,method), COLORS(method))
end
title('Time - rygcbmk')
