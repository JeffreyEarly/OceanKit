% BSplineBasisExample
%
% Visualize a B-spline basis and its first derivative basis.

K = 3;
t = (0:10)';
tKnot = [repmat(t(1), K-1, 1); t; repmat(t(end), K-1, 1)];
nSplines = numel(t) + K - 2;

basisSpline = BSpline(S=K-1, knotPoints=tKnot, xi=zeros(nSplines, 1));
tq = linspace(min(t), max(t), 1000)';

figure(Position=[100 100 860 560])
tiledlayout(2, 1, TileSpacing="compact")

nexttile
hold on
for iSpline = 1:nSplines
    coefficients = zeros(nSplines, 1);
    coefficients(iSpline) = 1;
    basisSpline.xi = coefficients;
    plot(tq, basisSpline(tq), LineWidth=1.5)
end
grid on
title("B-spline basis functions")

nexttile
hold on
for iSpline = 1:nSplines
    coefficients = zeros(nSplines, 1);
    coefficients(iSpline) = 1;
    basisSpline.xi = coefficients;
    plot(tq, basisSpline.valueAtPoints(tq, D=1), LineWidth=1.5)
end
grid on
xlabel("t")
title("First derivative basis functions")
