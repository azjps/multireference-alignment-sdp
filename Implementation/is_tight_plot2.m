
FONTSIZE = 14;

DEBUG = true;           % print status
L = 3;                  % number of rotations
SigmaVal = 1*[0:0.2:3];              % noise level
Nval = 1*[12:2:12];  % number of input signals
RUNS = 200;               % number of iterations to run over

METHODS = 1;

NvalNum  = length(Nval);
SigmaNum = length(SigmaVal);
resultsE = zeros(SigmaNum, METHODS);      % estimator error results
times    = zeros(SigmaNum, METHODS);      % runtime

% returns (positive) complex phase of z, normalized to be between [0,1)
add_one_if_neg = @(z) z + (z < 0);
pos_phase = @(z) add_one_if_neg(atan2(imag(z),real(z))/(2*pi));


for SigmaCounter = 1:SigmaNum
for run = 1:RUNS
    N = Nval(1);             
    sigma = SigmaVal(SigmaCounter); 
    
    x = random('normal',0,1,L,1);  % each component with size O(1);
    x = sqrt(L) * x / norm(x);
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
        shift = 0; % random('unid',L);
        trueShifts(k) = shift;
        Y(:,k) = circshift(x,shift) + sigma*random('normal',0,1,L,1);
        
        if DEBUG == true
             fprintf('Rotation %2d: ', k);
             disp(Y(:,k)');
        end
    end


    % run method and time
    tic
    xest = is_tight_SDPugJustPhase_pos(Y, sigma, true);
    runtime = toc;
    times(SigmaCounter,1) = times(SigmaCounter,1) + runtime;
    
    resultsE(SigmaCounter,1) = resultsE(SigmaCounter,1) + xest;

   
    if mod(run, 5) == 0
        fprintf('Run: %d, N: %d\n', run, N);
    end
end

clf 

% average and plot
resultsE(SigmaCounter,:) = resultsE(SigmaCounter,:) / RUNS;
times(SigmaCounter,:) = times(SigmaCounter,:) / RUNS;

%figure(73);
COLORS = 'rgcbmky';

subplot(1,1,1)
hold on

for method = 1:METHODS
    plot(SigmaVal(1:SigmaCounter), resultsE(1:SigmaCounter,method), COLORS(method), 'LineWidth', 2)
end
title('Tightness of Align-SDP+, N = 4, L = 3', 'FontSize', FONTSIZE)
xlabel('Sigma', 'FontSize', FONTSIZE)
ylabel('Ratio of Align-SDP+ Instances which are Tight', 'FontSize', FONTSIZE)

% subplot(1,2,2)
% hold on
% for method = 1:METHODS
%     semilogy(Nval(1:SigmaCounter), times(1:SigmaCounter,method), COLORS(method))
% end
% %title('Time - rygcbmk')    
% title('Runtime')
% xlabel('Sigma')
% ylabel('Percentage of SDP1 instances which are not tight')

drawnow

end
   