clf

L = 1000;

%# load image
load spine
img = X;
figure(3)
subplot(121), imshow(img, map), axis on

%# convert pixel coordinates from cartesian to polar
[r c] = size(img)
zcenter = r/2 + c/2*1i
rs = 1:1:round(min(r,c)/2-20)
[pim,zoff] = PolarTransform(img,zcenter,rs,2*pi/L);

pim2 = pim(:);
zoff2 = zoff(:);
zs = size(zoff);

subplot(122), warp(imag(zoff), real(zoff), zeros(size(zoff)), pim, map), axis ij square
view(2)


% figure(2)
% hold on
% plot(zoff,'o');



% [X Y] = meshgrid(1:c,1:r);
% [theta rho] = cart2pol(X, Y);

% %# show pixel locations (subsample to get less dense points)
% XX = X(1:8:end,1:4:end);
% YY = Y(1:8:end,1:4:end);
% tt = theta(1:8:end,1:4:end);
% rr = rho(1:8:end,1:4:end);
% subplot(121), scatter(XX(:),YY(:),3,'filled')
% subplot(122), scatter(tt(:),rr(:),3,'filled')
% 
% %# show images
% figure
% subplot(121), imshow(img, map), axis on
% subplot(122), warp(theta, rho, zeros(size(theta)), img, map)
% view(2), axis square

