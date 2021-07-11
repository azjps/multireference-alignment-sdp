
DEBUG = true;           % print status
L = 12;                  % number of rotations
sigma = 1.5;              % noise level
Nval = 1*[10:10:100];  % number of input signals
RUNS = 200;               % number of iterations to run over

estimator_list = {                               ...
    @pickfirst,                                  ... % Pick first
    @InvarianceMethod,                           ... % Invariance method
    @angsynchshiftspec,                          ... % Do spectral angular synchornization for the shifts
    @UGJustPhase_kmeans,                        ... % k-means
}; 

estimator_names = {'First','Invariants','Angular Sync','SDP-'};

ShiftRecovery(DEBUG, L, sigma, Nval, RUNS, estimator_list, estimator_names);

%    @SDPugJustPhase_pos,                         ... % SDP MLE Synch recovering just shifts  


%    @InvarianceMethod,                           ... % Invariance method
%    @angsynchshiftspec,                          ... % Do spectral angular synchornization for the shifts
%    @UGJustPhase_kmeans,                        ... % k-means

%    @UGJustPhase_opt_short,                      ... % SDP-less method for MLE, round best signal to rotation 

%    @UGJustPhase_opt_short,                     ... % SDP-less method for MLE, round first signal to rotation 1                   
%    @UGJustPhase,                               ... % SDP-less method for MLE
%    @SDPugJustPhase                             ... % SDP MLE Synch recovering just shifts  
%    @SPECTRALugJustPhase,                       ... % Spectral MLE Synch recovering just shifts      
%    @SDPforMLEshiftandsignalMoreConstraints,    ... % SDP MLE Synch recovering signal, using more constraints 
%    @SDPforMLEshiftandsignal,                   ... % SDP MLE Synch recovering signal too
%    @SDPforMLEshiftandsignalNonIndicators,      ... % SDP MLE Synch recovering signal too, non-indicator
