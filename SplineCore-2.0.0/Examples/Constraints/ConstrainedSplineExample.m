% ConstrainedSplineExample
%
% Choose one knot pattern and see how knot placement changes a simple
% one-dimensional constrained fit. Set `knotType` to one of:
%
% 1. automatic knots
% 2. transition knots
% 3. wider transition knots
% 4. single repeated midpoint knot
% 5. repeated midpoint knots
% 6. repeated transition knots

knotType = 6;

switch knotType
    case 1
        knotPoints = [];
        useTransitionConstraints = false;
        summary = "Automatic unconstrained fit.";
    case 2
        knotPoints = [0; 4; 4.5; 5.5; 6; 10];
        useTransitionConstraints = true;
        summary = "Transition knots with local stationary and zero-acceleration constraints.";
    case 3
        knotPoints = [0; 3.5; 4.5; 5.5; 6.5; 10];
        useTransitionConstraints = true;
        summary = "Wider transition region with the same local constraints.";
    case 4
        knotPoints = [0; 5; 5; 10];
        useTransitionConstraints = true;
        summary = "Single repeated midpoint knot with the same local constraints.";
    case 5
        knotPoints = [0; 4; 4; 6; 6; 10];
        useTransitionConstraints = true;
        summary = "Repeated midpoint knots that allow sharper regime changes.";
    case 6
        knotPoints = [0; 4.5; 4.5; 5.5; 5.5; 10];
        useTransitionConstraints = true;
        summary = "Repeated transition knots centered on the observed jump.";
    otherwise
        error("ConstrainedSplineExample:UnknownCase", "knotType must be an integer from 1 to 6.");
end

t = (0:10)';
x = [0; 0; 0; 0; 0; 0; 2; 4; 6; 8; 10];
tq = linspace(t(1), t(end), 1000)';

S = 2;
if useTransitionConstraints
    constraints = PointConstraint.equal([2.5; 2.5; 7.5], D=[1; 2; 2], value=0);
else
    constraints = PointConstraint.empty(0,1);
end

fit = ConstrainedSpline(t, x, S=S, knotPoints=knotPoints, constraints=constraints);

figure(Position=[100 100 900 700])
tiledlayout(3, 1, TileSpacing="compact")

panelLabels = ["Position", "1st derivative", "2nd derivative"];
for derivativeOrder = 0:S
    nexttile
    plot(tq, fit.valueAtPoints(tq, D=derivativeOrder), LineWidth=2), hold on

    if derivativeOrder == 0
        scatter(t, x, 45, "filled")
        ylim([-1 11])
    else
        addConstraintLines(constraints, derivativeOrder)
    end

    addKnotLines(knotPoints)
    xlim([t(1) t(end)])
    grid on
    ylabel(panelLabels(derivativeOrder + 1))

    if derivativeOrder < S
        set(gca, XTickLabel=[])
    else
        xlabel("Time")
    end
end

sgtitle(sprintf("Case %d: %s", knotType, summary))

function addKnotLines(knotPoints)
if isempty(knotPoints)
    return
end

for knotValue = unique(knotPoints(:)).'
    xline(knotValue, ":", Color=[0.2 0.6 0.2], LineWidth=1);
end
end

function addConstraintLines(constraints, derivativeOrder)
if isempty(constraints)
    return
end

matchingConstraints = constraints(arrayfun(@(c) isequal(c.D, derivativeOrder), constraints));
if isempty(matchingConstraints)
    return
end

for point = unique(vertcat(matchingConstraints.points)).'
    xline(point, "--", Color=[0.85 0.2 0.2], LineWidth=1);
end
end
