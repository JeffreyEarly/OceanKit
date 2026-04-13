%% Tutorial Metadata
% Title: Robust Fitting with Outliers
% Slug: robust-fitting-with-outliers
% Description: Replace a normal noise model with a Student-t model when a few observations are badly contaminated.
% NavOrder: 3

%% Switch the noise model from normal to Student-t
% The workflow is the same as in the previous
% [`ConstrainedSpline`](../classes/constrainedspline) tutorial, but now a
% few observations are badly contaminated. The key change is the
% [`distribution`](../classes/constrainedspline/distribution) object: use
% a [`StudentTDistribution`](https://en.wikipedia.org/wiki/Student%27s_t-distribution)
% when you want the fit to be less sensitive to outliers than a
% [normal model](https://en.wikipedia.org/wiki/Normal_distribution).
%
% $$
% y_i = f(t_i) + \epsilon_i, \qquad \epsilon_i \sim t_{\nu}(0,\sigma).
% $$

rng(7)
t = linspace(0, 1, 60)';
xTrue = exp(-3*t).*sin(4*pi*t);
xObs = xTrue + 0.08*randn(size(t));
outlierIndex = [10 22 37 51];
xObs(outlierIndex) = xObs(outlierIndex) + [0.75; -0.55; 0.65; -0.70];
tq = linspace(t(1), t(end), 400)';

normalModel = NormalDistribution(sigma=0.08);
studentTModel = StudentTDistribution(sigma=0.08, nu=3);

normalFit = ConstrainedSpline.fromData(t, xObs, S=3, splineDOF=12, distribution=normalModel);
robustFit = ConstrainedSpline.fromData(t, xObs, S=3, splineDOF=12, distribution=studentTModel);

% Plot the normal and Student-t fits on the same contaminated data.
figure(Position=[100 100 900 320])
plot(tq, exp(-3*tq).*sin(4*pi*tq), "k--", LineWidth=1.5), hold on
plot(tq, normalFit(tq), LineWidth=2)
plot(tq, robustFit(tq), LineWidth=2)
scatter(t, xObs, 28, "filled", MarkerFaceAlpha=0.65)
scatter(t(outlierIndex), xObs(outlierIndex), 72, "o", LineWidth=1.5)
xlabel("t")
ylabel("x(t)")
legend("Underlying signal", "Normal model", "Student-t model", "Observations", "Tagged outliers", Location="northeast")
grid on
if exist("tutorialFigureCapture", "var") && isa(tutorialFigureCapture, "function_handle"), tutorialFigureCapture("robust-fit-comparison", Caption="Changing the distribution from normal to Student-t makes the fit less sensitive to a few large outliers."); end

%% Inspect the final robust weights
% Internally the robust fit is solved by an iteratively reweighted scheme,
% but the tutorial-level idea is simply that the Student-t model downweights
% observations that look much less plausible under a normal model. For the
% implementation details, see
% [`distribution`](../classes/constrainedspline/distribution)
% and the developer topic
% [`tensorModelSolution`](../classes/constrainedspline/tensormodelsolution).

% Plot the final weights assigned by the robust fit.
figure(Position=[100 100 780 300])
stem(t, robustFit.W, "filled", LineWidth=1.2), hold on
scatter(t(outlierIndex), robustFit.W(outlierIndex), 70, "o", LineWidth=1.5)
xlabel("t")
ylabel("Final weight")
grid on
if exist("tutorialFigureCapture", "var") && isa(tutorialFigureCapture, "function_handle"), tutorialFigureCapture("robust-fit-weights", Caption="The final weights show which observations the Student-t fit decided to trust less."); end
