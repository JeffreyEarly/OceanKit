% MonotonicSplineExample
%
% Compare unconstrained and globally monotone fits on noisy data.

rng(1)
distribution = StudentTDistribution(sigma=0.2, nu=3.0);

t = linspace(0, 2, 11)';
f = @(x) tanh(x - 1);
x = f(t) + distribution.rand(size(t));

S = 2;
knotPoints = [min(t); min(t); 0.5; 1; 1.5; max(t); max(t)];
tq = linspace(min(t), max(t), 400)';

unconstrainedSpline = ConstrainedSpline.fromData(t, x, S=S, knotPoints=knotPoints, distribution=distribution);

monotoneSpline = ConstrainedSpline.fromData(t, x, S=S, knotPoints=knotPoints, distribution=distribution, constraints=GlobalConstraint.monotonicIncreasing());

figure(Position=[100 100 860 560])
tiledlayout(2, 1, TileSpacing="compact")

nexttile
plot(tq, unconstrainedSpline(tq), LineWidth=2), hold on
plot(tq, monotoneSpline(tq), LineWidth=2)
scatter(t, x, 45, "filled")
grid on
ylabel("Value")
title("Monotonic spline fit")
legend("Unconstrained", "Monotone increasing", "Samples",  Location="southoutside")

nexttile
plot(tq, unconstrainedSpline.valueAtPoints(tq, D=1), LineWidth=2), hold on
plot(tq, monotoneSpline.valueAtPoints(tq, D=1), LineWidth=2)
yline(0, "k--")
grid on
xlabel("Time")
ylabel("First derivative")
title("The global monotonicity constraint keeps the derivative nonnegative")
