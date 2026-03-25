% SplineDifferentiationExample
%
% Compare direct derivative evaluation with the derivative spline returned
% by diff.

f = @(x) sin(4*pi*x);
df = @(x) 4*pi*cos(4*pi*x);
d3f = @(x) -(4*pi)^3*cos(4*pi*x);

N = 11;
S = 3;
x = linspace(-1, 1, N)';
xDense = linspace(-1, 1, 400)';

splineFit = InterpolatingSpline(x, f(x), S=S);
derivativeSpline = diff(splineFit, 1);

figure(Position=[100 100 980 560])
tiledlayout(2, 2, TileSpacing="compact")

nexttile
plot(xDense, f(xDense), "k--", LineWidth=1.5), hold on
plot(xDense, splineFit(xDense), LineWidth=2)
scatter(x, f(x), 36, "filled")
grid on
title("Interpolated signal")
legend("Truth", "Spline", "Samples", Location="southoutside")

nexttile
plot(xDense, df(xDense), "k--", LineWidth=1.5), hold on
plot(xDense, splineFit.valueAtPoints(xDense, D=1), LineWidth=2)
grid on
title("Direct first derivative evaluation")
legend("Truth", "Spline derivative", Location="southoutside")

nexttile
plot(xDense, df(xDense), "k--", LineWidth=1.5), hold on
plot(xDense, derivativeSpline(xDense), LineWidth=2)
grid on
title("Derivative spline from diff")
legend("Truth", "diff(spline)", Location="southoutside")

nexttile
plot(xDense, d3f(xDense), "k--", LineWidth=1.5), hold on
plot(xDense, splineFit.valueAtPoints(xDense, D=3), LineWidth=2)
plot(xDense, derivativeSpline.valueAtPoints(xDense, D=2), ":", LineWidth=2)
grid on
title("Higher derivatives remain consistent")
legend("Truth", "Direct D = 3", "diff(spline), D = 2", Location="southoutside")
