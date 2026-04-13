% RungeInterpolationExample
%
% Compare uniform and Chebyshev sampling when interpolating the Runge
% function.

f = @(x) 1./(25*x.^2 + 1);
N = 11;
S = 3;
xDense = linspace(-1, 1, 400)';

xUniform = linspace(-1, 1, N)';
uniformSpline = InterpolatingSpline.fromGriddedValues(xUniform, f(xUniform), S=S);
uniformError = max(abs(f(xDense) - uniformSpline(xDense)));

xChebyshev = cos((0:N-1)'*pi/(N-1));
chebyshevSpline = InterpolatingSpline.fromGriddedValues(xChebyshev, f(xChebyshev), S=S);
chebyshevError = max(abs(f(xDense) - chebyshevSpline(xDense)));

xReference = [-1; -0.33; 0.33; 1];
yReference = [0; 1; 2; 0];
referenceSpline = InterpolatingSpline.fromGriddedValues(xReference, yReference, S=3);
referenceInterpolant = griddedInterpolant(xReference, yReference, "spline");

figure(Position=[100 100 980 620])
tiledlayout(2, 2, TileSpacing="compact")

nexttile
plot(xDense, f(xDense), "k--", LineWidth=1.5), hold on
plot(xDense, uniformSpline(xDense), LineWidth=2)
scatter(xUniform, f(xUniform), 36, "filled")
title(sprintf("Uniform grid, max error %.2g", uniformError))
grid on

nexttile
plot(xDense, f(xDense), "k--", LineWidth=1.5), hold on
plot(xDense, chebyshevSpline(xDense), LineWidth=2)
scatter(xChebyshev, f(xChebyshev), 36, "filled")
title(sprintf("Chebyshev grid, max error %.2g", chebyshevError))
grid on

nexttile([1 2])
scatter(xReference, yReference, 45, "filled"), hold on
plot(xDense, referenceSpline(xDense), "b", LineWidth=2.5)
plot(xDense, referenceInterpolant(xDense), "k--", LineWidth=1.5)
grid on
xlabel("x")
ylabel("y")
title("InterpolatingSpline agrees with MATLAB spline interpolation in 1D")
legend("Samples", "InterpolatingSpline", "griddedInterpolant", Location="southoutside")
