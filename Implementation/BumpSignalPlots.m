PLOTX = 3;
PLOTY = 1;
FONTSIZE = 16;
                                             
x = [0:0.01:1];
y_clean = zeros(101,1);
y_clean([10:35]) = ones(26,1);

y_shifted = zeros(101,1);
y_shifted([60:85]) = ones(26,1);

y_shifted2 = zeros(101,1);
y_shifted2([30:55]) = ones(26,1);

y_shifted_noisy = y_shifted + 0.25*randn(101,1);
y_shifted2_noisy = y_shifted2 + 0.25*randn(101,1);

y_clean_noisy = y_clean + 0.25*randn(101,1);

y_clean_very_noisy = y_clean + randn(101,1);
y_shifted_very_noisy = y_shifted + randn(101,1);
y_shifted2_very_noisy = y_shifted2 + randn(101,1);

figure(2);
subplot(PLOTX,PLOTY,1);
set(gca, 'FontSize', FONTSIZE)
plot(x,y_clean,'LineWidth',4,'Color','blue')
axis([0 1 -1 1.5])
title('Clean Template','Fontsize', FONTSIZE);

subplot(PLOTX,PLOTY,2);
set(gca, 'FontSize', FONTSIZE)
plot(x,y_shifted,'LineWidth',4,'Color','blue')
axis([0 1 -1 1.5])
title('Shifted Template','Fontsize', FONTSIZE);

subplot(PLOTX,PLOTY,3);
set(gca, 'FontSize', FONTSIZE)
plot(x,y_shifted2,'LineWidth',4,'Color','blue')
axis([0 1 -1 1.5])
title('Shifted Template','Fontsize', FONTSIZE);

figure(3);
subplot(PLOTX,PLOTY,1);
set(gca, 'FontSize', FONTSIZE)
plot(x,y_clean,'LineWidth',4,'Color','blue')
axis([0 1 -1 1.5])
title('Clean Template','Fontsize', FONTSIZE);

subplot(PLOTX,PLOTY,2);
set(gca, 'FontSize', FONTSIZE)
plot(x,y_shifted_noisy,'LineWidth',4,'Color','green')
axis([0 1 -1 1.5])
title('Noisy Reference','Fontsize', FONTSIZE);

subplot(PLOTX,PLOTY,3);
set(gca, 'FontSize', FONTSIZE)
plot(x,y_shifted2_noisy,'LineWidth',4,'Color','green')
axis([0 1 -1 1.5])
title('Noisy Reference','Fontsize', FONTSIZE);

figure(4);
subplot(PLOTX,PLOTY,1);
set(gca, 'FontSize', FONTSIZE)
plot(x,y_clean_noisy,'LineWidth',4,'Color','green')
axis([0 1 -2 2.5])
title('Noisy Reference','Fontsize', FONTSIZE);

subplot(PLOTX,PLOTY,2);
set(gca, 'FontSize', FONTSIZE)
plot(x,y_shifted_noisy,'LineWidth',4,'Color','green')
axis([0 1 -2 2.5])
title('Noisy Reference','Fontsize', FONTSIZE);

subplot(PLOTX,PLOTY,3);
set(gca, 'FontSize', FONTSIZE)
plot(x,y_shifted2_noisy,'LineWidth',4,'Color','green')
axis([0 1 -2 2.5])
title('Noisy Reference','Fontsize', FONTSIZE);

figure(5);
subplot(PLOTX,PLOTY,1);
set(gca, 'FontSize', FONTSIZE)
plot(x,y_clean_very_noisy,'LineWidth',4,'Color','red')
axis([0 1 -2 2.5])
title('Very Noisy Reference','Fontsize', FONTSIZE);

subplot(PLOTX,PLOTY,2);
set(gca, 'FontSize', FONTSIZE)
plot(x,y_shifted_very_noisy,'LineWidth',4,'Color','red')
axis([0 1 -2 2.5])
title('Very Noisy Reference','Fontsize', FONTSIZE);

subplot(PLOTX,PLOTY,3);
set(gca, 'FontSize', FONTSIZE)
plot(x,y_shifted2_very_noisy,'LineWidth',4,'Color','red')
axis([0 1 -2 2.5])
title('Very Noisy Reference','Fontsize', FONTSIZE);
