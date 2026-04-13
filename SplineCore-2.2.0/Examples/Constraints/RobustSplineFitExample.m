% RobustSplineFitExample
%
% Compare ordinary least squares and a robust Student-t fit on noisy data.

rng(2)
distribution = StudentTDistribution(sigma=1, nu=3.0);

t = linspace(0, 2, 11)';
f = @(x) x.^2 + 2*x + 1;
x = f(t) + distribution.rand(size(t));
x([4 9]) = x([4 9]) + [2.5; -2.0];

S = 2;
knotPoints = [min(t); min(t); 0.75; 1.25; max(t); max(t)];
tq = linspace(min(t), max(t), 400)';

leastSquaresSpline = ConstrainedSpline.fromData(t, x, S=S, knotPoints=knotPoints);
robustSpline = ConstrainedSpline.fromData(t, x, S=S, knotPoints=knotPoints, distribution=distribution);

figure(Position=[100 100 860 520])
tiledlayout(2, 1, TileSpacing="compact")

nexttile
plot(tq, f(tq), "k--", LineWidth=1.5), hold on
plot(tq, leastSquaresSpline(tq), LineWidth=2)
plot(tq, robustSpline(tq), LineWidth=2)
scatter(t, x, 45, "filled")
grid on
ylabel("Value")
title("Robust spline fitting with outliers")
legend("Truth", "Least squares", "Student-t robust fit", "Samples",  Location="southoutside")

nexttile
stem(t, robustSpline.W, "filled", LineWidth=1.2), hold on
scatter(t([4 9]), robustSpline.W([4 9]), 70, "o", LineWidth=1.5)
grid on
xlabel("Time")
ylabel("Final IRLS weight")
title("Outliers receive lower weights in the robust fit")
