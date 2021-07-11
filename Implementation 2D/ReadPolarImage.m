
function ReadPolarImage

clf

L = 80;        % number of rotations to draw from
R = 80;        % number of radial discretizations
Nval = 40:1:40;         % number of pictures
RUNS = 1;       % number of iterations to average over
sigma = 4;     % level of noise
DEBUG = false;

%# load image
[img2, map] = rgb2ind(imread('peace_sign.jpg'), 128);
% hack for peace sign
[r c] = size(img2);
img = 4*ones(800, 800);
img(10:10+r-1,100:100+c-1) = img2;

%load spine
%img = X
figure(3)
subplot(331), imshow(img, map), axis off square

%# convert pixel coordinates from cartesian to polar
[r c] = size(img);
zcenter = r/2 + c/2*1i;
zradius = round(min(r,c)/2-1);
rs = 1:round(zradius/R):zradius;
R = size(rs);
R = R(2);
[pim,zoff] = PolarTransform(img,zcenter,rs,2*pi/L);
% pim is pixel intensity matrix, where pixel i,j
% at jth radial line, ith radial discretization
% so size(pim) = R x L
% so circshift by row to rotate in T
% is located at 

estimator_list = {                              ...
    @BestRotationWithFirst,                     ... 
}; 
%    @SDP_UG,                                    ... % SDP MLE Synch recovering just shifts  
%    @SDP_UG_pos,                                ... % SDP MLE Synch recovering just shifts  
METHODS = length(estimator_list);


NvalNum  = length(Nval);
resultsE = zeros(NvalNum, METHODS);      % estimator error results
times    = zeros(NvalNum, METHODS);      % runtime

for Ncounter = 1:NvalNum
    for run = 1:RUNS
        
        N = Nval(Ncounter);   
        fprintf('Run #: %d, Number of images: %d\n', run, N);
        
        Y = zeros(R,L,N);
        for i = 1:N
           Y(:,:,i) = RotateImage(pim, random('unid', L), sigma);            
        end
        
        for method = 1:METHODS
            % run method and time
            tic
            Yest = estimator_list{method}(Y, sigma, DEBUG);
            runtime = toc;
            times(Ncounter,method) = times(Ncounter,method) + runtime;

            error = ShiftLessDist(pim, Yest);
            if DEBUG == true
                 fprintf('Estimate %d with error %2.3f, time %2.3f: \n', ...
                     [method; error; runtime]);
                 disp(Yest');
            end
            resultsE(Ncounter,method) = resultsE(Ncounter,method) + error;
            
            if Ncounter == 1 && run == 1 && method == 1                
         
                subplot(332), warp(imag(zoff), real(zoff), zeros(size(zoff)), RotateImage(pim, 0, 0), map), axis off square
                view(2)

                subplot(334), warp(imag(zoff), real(zoff), zeros(size(zoff)), RotateImage(pim, round(L/5), sigma), map), axis off square
                view(2)

                subplot(335), warp(imag(zoff), real(zoff), zeros(size(zoff)), RotateImage(Yest, 0, 0), map), axis off square
                view(2)
                drawnow 
            elseif Ncounter == 1 && run == 1 && method == 2
                subplot(336), warp(imag(zoff), real(zoff), zeros(size(zoff)), RotateImage(Yest, 0, 0), map), axis off square
                view(2)
                drawnow 
            end
        end
    end    
    
% average and plot
resultsE(Ncounter,:) = resultsE(Ncounter,:) / RUNS;
times(Ncounter,:) = times(Ncounter,:) / RUNS;

figure(4)
COLORS = 'rgcbmky';

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
%title('Time - rygcbmk')    
title('rygcbmk')

drawnow    
    
end

end