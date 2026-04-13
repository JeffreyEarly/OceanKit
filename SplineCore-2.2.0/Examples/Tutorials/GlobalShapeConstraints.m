%% Tutorial Metadata
% Title: Global Shape Constraints
% Slug: global-shape-constraints
% Description: Enforce positivity and monotonicity over an entire one-dimensional domain.
% NavOrder: 5

%% Enforce positivity and monotonicity across the whole domain
% [`GlobalConstraint`](../classes/constraints/globalconstraint) objects do
% not act at a few selected points. Instead they describe the shape of a
% [`ConstrainedSpline`](../classes/constrainedspline) fit everywhere on
% the domain:
%
% $$
% f(t)\ge 0,\qquad f'(t)\ge 0 \quad \text{for all } t.
% $$

rng(4)
t = linspace(0, 1, 45)';
xTrue = 0.15 + 0.85*(1 - exp(-4*t));
xObs = xTrue + 0.05*randn(size(t));
tq = linspace(t(1), t(end), 400)';

noiseModel = NormalDistribution(sigma=0.05);
freeFit = ConstrainedSpline.fromData(t, xObs, S=3, splineDOF=12, distribution=noiseModel);

shapeConstraints = [ ...
    GlobalConstraint.positive()
    GlobalConstraint.monotonicIncreasing()];

shapeFit = ConstrainedSpline.fromData(t, xObs, S=3, splineDOF=12, distribution=noiseModel, constraints=shapeConstraints);

% Plot the unconstrained and globally constrained fits together.
figure(Position=[100 100 820 320])
plot(tq, 0.15 + 0.85*(1 - exp(-4*tq)), "k--", LineWidth=1.5), hold on
plot(tq, freeFit(tq), LineWidth=2)
plot(tq, shapeFit(tq), LineWidth=2)
scatter(t, xObs, 28, "filled", MarkerFaceAlpha=0.65)
xlabel("t")
ylabel("x(t)")
legend("Underlying signal", "Unconstrained fit", "Positive monotone fit", "Observations", Location="southoutside")
grid on
if exist("tutorialFigureCapture", "var") && isa(tutorialFigureCapture, "function_handle"), tutorialFigureCapture("positive-monotone-fit", Caption="GlobalConstraint objects can enforce positivity and monotonic increase over the full domain in one fit."); end

%% Check the constrained derivative
% The first derivative is the clearest way to see the monotonicity
% condition. Once the global monotone-increasing constraint is active, the
% fitted derivative stays nonnegative throughout the interval.

dFree = freeFit.valueAtPoints(tq, D=1);
dShape = shapeFit.valueAtPoints(tq, D=1);

% Plot the first derivative to check the monotonicity condition directly.
figure(Position=[100 100 780 300])
plot(tq, dFree, LineWidth=1.6), hold on
plot(tq, dShape, LineWidth=2)
yline(0, "k--")
xlabel("t")
ylabel("dx/dt")
legend("Unconstrained", "Positive monotone fit", Location="southoutside")
grid on
if exist("tutorialFigureCapture", "var") && isa(tutorialFigureCapture, "function_handle"), tutorialFigureCapture("positive-monotone-derivative", Caption="The constrained fit keeps the first derivative nonnegative across the whole domain."); end
