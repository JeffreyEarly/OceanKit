% BSplineFiniteDifferenceComparison
%
% Compare a spline-based second-derivative operator with a standard
% centered finite-difference approximation.

K = 3;
t = linspace(0, 20, 21)';
tKnot = BSpline.knotPointsForDataPoints(t, S=K-1);

B = BSpline.matrixForDataPoints(t, knotPoints=tKnot, S=K-1, D=2);
X = B(:,:,1);
A = B(:,:,3);
splineSecondDerivativeOperator = A / X;

dt = t(2) - t(1);
finiteDifferenceOperator = zeros(numel(t));
for iPoint = 2:numel(t)-1
    finiteDifferenceOperator(iPoint, iPoint-1:iPoint+1) = [1 -2 1] / dt^2;
end
finiteDifferenceOperator(1, 1:3) = [1 -2 1] / dt^2;
finiteDifferenceOperator(end, end-2:end) = [1 -2 1] / dt^2;

f = @(x) sin(0.35*x);
d2f = @(x) -(0.35^2)*sin(0.35*x);
values = f(t);

splineAcceleration = splineSecondDerivativeOperator * values;
finiteDifferenceAcceleration = finiteDifferenceOperator * values;

figure(Position=[100 100 980 560])
tiledlayout(2, 2, TileSpacing="compact")

nexttile
imagesc(splineSecondDerivativeOperator)
axis image
colorbar
title("Spline second-derivative operator")

nexttile
imagesc(finiteDifferenceOperator)
axis image
colorbar
title("Finite-difference operator")

nexttile([1 2])
plot(t, d2f(t), "k--", LineWidth=1.5), hold on
plot(t, splineAcceleration, "o-", LineWidth=1.5)
plot(t, finiteDifferenceAcceleration, "s-", LineWidth=1.5)
grid on
xlabel("t")
ylabel("Second derivative")
legend("Truth", "Spline operator", "Finite difference", Location="southoutside")
