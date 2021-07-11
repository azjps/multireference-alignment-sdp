function ShiftLessRecovery

clf

L = 3;                  % number of rotations
sigma = 0;              % noise level
Nval = 1*[5:5:35];      % number of input signals
NMax0 = 35;
Nmethods = 5;           % number of different techniques to try
Nruns = 20;             % number of iterations to run over

DEBUG = true;           % print status

NvalNum = length(Nval);
resultsE = zeros(NvalNum,Nmethods);
times = zeros(NvalNum,Nmethods);

for run = 1:Nruns
for Ncounter = 1:NvalNum
    N = Nval(Ncounter);             
    
    x = random('normal',0,1,L,1);  % each component with size O(1);

    % measurements
    Y = zeros(L,N);

    trueShifts = zeros(N,1);
    
    for k = 1:N
        shift = random('unid',L);
        trueShifts(k) = shift;
        Y(:,k) = circshift(x,shift) + sigma*random('normal',0,1,L,1);
        
        if DEBUG == true
             fprintf('Rotation %2d: %f', [k; Y(:k)])
        end
    end


    %ON just pick the first one!
        method = 1;
        
        tic
        xest = pickfirst(Y,sigma);
        error = distshiftless(x,xest);
        resultsE(Ncounter,method) = resultsE(Ncounter,method) + error;
        times(Ncounter,method) = times(Ncounter,method) + toc;
        
    %ON SDP MLE Synch recovering just shifts
        method = 2;
        
        tic
        xest = SDPugJustPhase(Y, sigma);
        %xest = pickfirst(Y,sigma);
        error = distshiftless(x,xest);
        resultsE(Ncounter,method) = resultsE(Ncounter,method) + error;
        times(Ncounter,method) = times(Ncounter,method) + toc;
        
    %ON do spectral angular synchornization for the shifts
        method = 3;
        
        tic
        xest = angsynchshiftspec(Y,sigma);
        %xest = pickfirst(Y,sigma);
        error = distshiftless(x,xest);
        resultsE(Ncounter,method) = resultsE(Ncounter,method) + error;
        times(Ncounter,method) = times(Ncounter,method) + toc;
   
    %OFF SDP MLE Synch recovering signal too
        method = 4;
        
        tic
        xest = SDPforMLEshiftandsignal(Y,sigma);
        %xest = pickfirst(Y,sigma);
        error = distshiftless(x,xest);
        resultsE(Ncounter,method) = resultsE(Ncounter,method) + error;
        times(Ncounter,method) = times(Ncounter,method) + toc;
        
        
    %OFF SDP MLE Synch recovering signal too, non-indicator
        method = 5;
        
        tic
        xest = SDPforMLEshiftandsignalNonIndicators(Y,sigma);
        %xest = pickfirst(Y,sigma);
        error = distshiftless(x,xest);
        resultsE(Ncounter,method) = resultsE(Ncounter,method) + error;
        times(Ncounter,method) = times(Ncounter,method) + toc;
    
   
fprintf('\n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n' )
disp(['Completed: N = ',int2str(N),' out of ',int2str(NMax0)])
disp(['On run %: ',int2str(floor(100*run/Nruns)),'% '])
fprintf('\n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n' )
        
end
fprintf('\n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n' )
disp(['Completed: ',int2str(floor(100*run/Nruns)),'% '])
fprintf('\n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n' )
end



resultsE = resultsE/Nruns
times = times/Nruns


figure(73)

subplot(1,2,1)
hold on
%method 1
    p1 = resultsE(:,1);
    plot(Nval,p1,'y')
%method 2
    p1 = resultsE(:,2);
    plot(Nval,p1,'g')
%method 3
    p1 = resultsE(:,3);
    plot(Nval,p1,'r')
%method 4
    p1 = resultsE(:,4);
    plot(Nval,p1,'k')
%method 5
    p1 = resultsE(:,5);
    plot(Nval,p1,'b')
title('error varying with N')


subplot(1,2,2)
hold on
%method 1
    p1 = times(:,1);
    plot(Nval,p1,'y')
%method 2
    p1 = times(:,2);
    plot(Nval,p1,'g')
%method 3
    p1 = times(:,3);
    plot(Nval,p1,'r')
%method 4
    p1 = times(:,4);
    plot(Nval,p1,'k')
%method 5
    p1 = times(:,5);
    plot(Nval,p1,'b')
title('time - 1y 2g 3r 4k 5b')


%disp('little test')
%Origx = x
%Corrx = PowerSpectrumCorrection(x,Y,sigma)



end