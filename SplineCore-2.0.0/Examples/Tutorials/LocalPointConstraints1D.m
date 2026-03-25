%% Tutorial Metadata
% Title: Local Point Constraints
% Slug: local-point-constraints
% Description: Apply value, slope, and curvature constraints at specific points in a one-dimensional spline fit.
% NavOrder: 4

%% Constrain the value, slope, or curvature at one point
% Local constraints act at specific coordinates. The derivative-order
% option `D` tells
% [`PointConstraint.equal`](../classes/constraints/pointconstraint/equal)
% whether you are constraining the value, the slope, or a higher
% derivative in a [`ConstrainedSpline`](../classes/constrainedspline) fit:
%
% $$
% f(t_c)=x_c,\qquad f'(t_c)=0,\qquad f''(t_c)=0.
% $$

rng(9)
t = linspace(0, 1, 32)';
xTrue = 0.2 + 0.45*sin(2*pi*t) + 0.18*cos(5*pi*t);
xObs = xTrue + 0.05*randn(size(t));
tq = linspace(t(1), t(end), 500)';

noiseModel = NormalDistribution(0.05);
freeFit = ConstrainedSpline(t, xObs, S=3, splineDOF=11, distribution=noiseModel);

valueConstraint = PointConstraint.equal(0.42, value=0.55);
flatConstraint = PointConstraint.equal(0.42, D=1, value=0);
curvatureConstraint = PointConstraint.equal(0.68, D=2, value=0);

valueFit = ConstrainedSpline(t, xObs, S=3, splineDOF=11, distribution=noiseModel, constraints=valueConstraint);
flatFit = ConstrainedSpline(t, xObs, S=3, splineDOF=11, distribution=noiseModel, constraints=flatConstraint);
curvatureFit = ConstrainedSpline(t, xObs, S=3, splineDOF=11, distribution=noiseModel, constraints=curvatureConstraint);

% Plot three single-constraint fits against the same unconstrained baseline.
caseFits = {valueFit, flatFit, curvatureFit};
caseTitles = ["Value target", "Zero slope", "Zero curvature"];
casePoints = [0.42 0.42 0.68];
caseTargets = [0.55 NaN NaN];

figure(Position=[100 100 980 760])
tiledlayout(2, 2, TileSpacing="compact")

for iCase = 1:3
    nexttile
    plot(tq, freeFit(tq), "--", LineWidth=1.5), hold on
    plot(tq, caseFits{iCase}(tq), LineWidth=2)
    scatter(t, xObs, 20, "filled", MarkerFaceAlpha=0.55)
    if isnan(caseTargets(iCase))
        xline(casePoints(iCase), ":", Color=[0.82 0.2 0.2], LineWidth=1.2)
    else
        scatter(casePoints(iCase), caseTargets(iCase), 70, "o", LineWidth=1.4, MarkerEdgeColor=[0.82 0.2 0.2])
    end
    title(caseTitles(iCase))
    xlabel("t")
    ylabel("x(t)")
    grid on
end

nexttile
plot(tq, freeFit(tq), "--", LineWidth=1.5), hold on
plot(tq, valueFit(tq), LineWidth=2)
plot(tq, flatFit(tq), LineWidth=2)
plot(tq, curvatureFit(tq), LineWidth=2)
scatter(t, xObs, 20, "filled", MarkerFaceAlpha=0.55)
title("Compare the constrained fits")
xlabel("t")
ylabel("x(t)")
grid on
legend("Unconstrained", "Value", "Slope", "Curvature", "Observations", Location="southoutside")
if exist("tutorialFigureCapture", "var") && isa(tutorialFigureCapture, "function_handle"), tutorialFigureCapture("single-constraint-cases", Caption="PointConstraint.equal can enforce a target value, a zero slope, or a zero curvature at selected points."); end

%% Combine several local conditions in one fit
% A constrained fit can enforce several pointwise conditions at once. The
% simplest pattern is to stack multiple
% [`PointConstraint`](../classes/constraints/pointconstraint) objects into
% the same `constraints` array.

combinedConstraints = [ ...
    PointConstraint.equal(0.42, value=0.55)
    PointConstraint.equal(0.42, D=1, value=0)];

combinedFit = ConstrainedSpline(t, xObs, S=3, splineDOF=11, distribution=noiseModel, constraints=combinedConstraints);
combinedSlope = combinedFit.valueAtPoints(tq, D=1);

% Plot the value constraint together with the resulting zero-slope condition.
figure(Position=[100 100 820 680])
tiledlayout(2, 1, TileSpacing="compact")

nexttile
plot(tq, freeFit(tq), "--", LineWidth=1.6), hold on
plot(tq, combinedFit(tq), LineWidth=2)
scatter(t, xObs, 24, "filled", MarkerFaceAlpha=0.65)
scatter(0.42, 0.55, 70, "o", LineWidth=1.4, MarkerEdgeColor=[0.82 0.2 0.2])
xline(0.42, ":", Color=[0.82 0.2 0.2], LineWidth=1.2)
xlabel("t")
ylabel("x(t)")
title("Value and zero-slope constraint at the same point")
grid on
legend("Unconstrained", "Constrained", "Observations", "Value target", Location="southoutside")

nexttile
plot(tq, combinedSlope, LineWidth=2), hold on
yline(0, "k--")
xline(0.42, ":", Color=[0.82 0.2 0.2], LineWidth=1.2)
xlabel("t")
ylabel("df/dt")
grid on
if exist("tutorialFigureCapture", "var") && isa(tutorialFigureCapture, "function_handle"), tutorialFigureCapture("combined-local-constraints", Caption="Several PointConstraint objects can be combined to enforce multiple local conditions in one fit."); end
