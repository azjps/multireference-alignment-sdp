
clf

L = 40;
Big = 100;
Diff = 5;
X = 2.5;

signal1 = zeros(Big, 1);
signal2 = zeros(Big, 1);
signal3 = zeros(Big, 1);
signal4 = zeros(Big, 1);

signal1(Diff  :  L+Diff-1) = normrnd(0, X, L, 1);
signal2(L+2*Diff:2*L+2*Diff-1) = signal1(Diff : L+Diff-1);
signal3(Diff  :  L+Diff-1) = signal1(Diff  :  L+Diff-1) + normrnd(0, X, L, 1);
signal4(L+2*Diff:2*L+2*Diff-1) = signal1(Diff : L+Diff-1) + normrnd(0, X, L, 1);

figure(2)
hold on
plot(1:Big, signal1, '--r', 'LineWidth', 2), axis off
plot(1:Big, signal2, '--b', 'LineWidth', 2), axis off

figure(3)
hold on
plot(1:Big, signal1, '--r', 'LineWidth', 2), axis off
plot(1:Big, signal4, '--b', 'LineWidth', 2), axis off

figure(4)
hold on
plot(1:Big, signal3, '--r', 'LineWidth', 2), axis off
plot(1:Big, signal4, '--b', 'LineWidth', 2), axis off
hold off