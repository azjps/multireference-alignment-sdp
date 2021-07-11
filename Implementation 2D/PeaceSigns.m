clf

L = 110;        % number of rotations to draw from
R = 100;        % number of radial discretizations
Nval = 40:1:40;         % number of pictures
RUNS = 1;       % number of iterations to average over
sigma = 4;     % level of noise
DEBUG = false;

plotX = 3;
plotY = 2;

%# load image
[img2, map] = rgb2ind(imread('peace_sign.jpg'), 128);
% hack for peace sign
[r c] = size(img2);
img = 4*ones(800, 800);
img(10:10+r-1,100:100+c-1) = img2;


[r c] = size(img);
zcenter = r/2 + c/2*1i;
zradius = round(min(r,c)/2-1);
rs = 1:round(zradius/R):zradius;
R = size(rs);
R = R(2);
[pim,zoff] = PolarTransform(img,zcenter,rs,2*pi/L);


rotations = [0; randi(L,8,1)];

figure(3)
for i = 1:plotX * plotY
    subplot(plotX,plotY,i), warp(imag(zoff), real(zoff), zeros(size(zoff)), RotateImage(pim, rotations(i), 0), map), axis off square
    view(2)
end

sigma = 1;
low_sigma_avg = zeros(R,L);
figure(4)
for i = 1 : plotX * plotY
    rpim = RotateImage(pim, rotations(i), sigma);
    % low_sigma_avg = low_sigma_avg + RotateImage(rpim, -rotations(i), 0);
    subplot(plotX,plotY,i), warp(imag(zoff), real(zoff), zeros(size(zoff)), rpim, map), axis off square
    view(2)
end

figure(5)
for i = 1 : plotX * plotY
    rpim = RotateImage(pim, 0, sigma);
    low_sigma_avg = low_sigma_avg + rpim;
    subplot(plotX,plotY,i), warp(imag(zoff), real(zoff), zeros(size(zoff)), rpim, map), axis off square
    view(2)
end
low_sigma_avg = low_sigma_avg / plotX / plotY;

sigma = 5;
high_sigma_avg = zeros(R,L);
figure(6)
for i = 1 : plotX * plotY
    rpim = RotateImage(pim, rotations(i), sigma);
    % low_sigma_avg = low_sigma_avg + RotateImage(rpim, -rotations(i), 0);
    subplot(plotX,plotY,i), warp(imag(zoff), real(zoff), zeros(size(zoff)), rpim, map), axis off square
    view(2)
end

figure(7)
for i = 1 : plotX * plotY
    rpim = RotateImage(pim, 0, sigma);
    high_sigma_avg = high_sigma_avg + rpim;
    subplot(plotX,plotY,i), warp(imag(zoff), real(zoff), zeros(size(zoff)), rpim, map), axis off square
    view(2)
end
high_sigma_avg = high_sigma_avg / plotX / plotY;

figure(8)
subplot(1,2,1), warp(imag(zoff), real(zoff), zeros(size(zoff)), low_sigma_avg, map), axis off square
view(2)
subplot(1,2,2), warp(imag(zoff), real(zoff), zeros(size(zoff)), high_sigma_avg, map), axis off square
view(2)