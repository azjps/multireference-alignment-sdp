function RelativeShifts

% configurable values
DEBUG = false;           % print status
L = 3;                  % number of rotations
sigma = 1;             % noise level
Nval = 1*[5:5:50];      % number of input signals
RUNS = 20;               % number of iterations to run over

estimator_list = {                              ... 
    @AngSync_Offsets,                            ... % SDP MLE Synch recovering just shifts
    @UG_Offsets                                  ...
};

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
    
    if DEBUG == true
        fprintf('\nRun number: %d, Number of sampled signals: %d\nUnderlying signal: ', ...
            [run; N]);
        disp(x');
    end
    
    % TODO: print out phases for each rotation of x
    
    % measurements
    Y = zeros(L,N);

    trueShifts = zeros(N,1);
    trueOffsets = zeros(N,N);

    for k = 1:N
        shift = random('unid',L);
        trueShifts(k) = shift;
        Y(:,k) = circshift(x,shift) + sigma*random('normal',0,1,L,1);
        
        if DEBUG == true
             fprintf('Rotation %2d: ', k);
             disp(Y(:,k)');
        end
    end
    
    for n1 = 1:N
        for n2 = 1:N
            trueOffsets(n1,n2) = mod(trueShifts(n1)-trueShifts(n2),L);
        end
    end
    
    disp(trueOffsets);
           
    for method = 1:METHODS
        % run method and time
        tic
        Zest = estimator_list{method}(x, Y, sigma, DEBUG);
        runtime = toc;
        resultsE(Ncounter,method) = resultsE(Ncounter,method) + nnz(Zest-trueOffsets)/(N^2);
        times(Ncounter,method) = times(Ncounter,method) + runtime;
        
        disp(Zest);
    end
    
end

clf 

% average and plot
resultsE(Ncounter,:) = resultsE(Ncounter,:) / RUNS;
times(Ncounter,:) = times(Ncounter,:) / RUNS;

%figure(73);
COLORS = 'rygcbmk';

subplot(1,2,1)
hold on

for method = 1:METHODS
    plot(Nval(1:Ncounter), resultsE(1:Ncounter,method), COLORS(method))
end
title('Error varying with N')

subplot(1,2,2)
hold on
for method = 1:METHODS
    plot(Nval(1:Ncounter), times(1:Ncounter,method), COLORS(method))
end
title('Time - rygcbmk')    

drawnow

end
   
