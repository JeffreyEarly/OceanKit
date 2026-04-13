% SplineEvaluationExample
%
% Compare three equivalent ways to evaluate a fitted spline and its
% derivatives: directly from the basis matrix, from piecewise-polynomial
% coefficients, and through the object-oriented interface.

f = @(x) sin(2*pi*x/10);

K = 3;
D = K - 1;
t = (0:10)';
tKnot = BSpline.knotPointsForDataPoints(t, S=K-1);
tqBasis = linspace(min(t), max(t), 1000)';
tqPP = linspace(min(t) - 1, max(t) + 1, 1000)';

B = BSpline.matrixForDataPoints(t, knotPoints=tKnot, S=K-1, D=D);
X = B(:,:,1);
coefficients = X\f(t);

Bq = BSpline.matrixForDataPoints(tqBasis, knotPoints=tKnot, S=K-1, D=D);
[C, tPP] = BSpline.ppCoefficientsFromSplineCoefficients(xi=coefficients, knotPoints=tKnot, S=K-1);
splineFit = InterpolatingSpline.fromGriddedValues(t, f(t), S=K-1);

titles = [
    "Evaluated with basis matrices"
    "Evaluated with PP coefficients"
    "Evaluated through InterpolatingSpline"
    ];

figure(Position=[100 100 980 780])
tiledlayout(K, 3, TileSpacing="compact")

for iColumn = 1:3
    for iDerivative = 0:D
        nexttile
        if iDerivative == 0
            scatter(t, f(t), 28, "filled"), hold on
        else
            hold on
        end

        switch iColumn
            case 1
                queryPoints = tqBasis;
                values = Bq(:,:,iDerivative + 1) * coefficients;
            case 2
                queryPoints = tqPP;
                values = BSpline.evaluateFromPPCoefficients(queryPoints=tqPP, C=C, tpp=tPP, D=iDerivative);
            case 3
                queryPoints = tqBasis;
                values = splineFit.valueAtPoints(tqBasis, D=iDerivative);
        end

        plot(queryPoints, values, LineWidth=1.75)
        grid on
        if iDerivative == 0
            title(titles(iColumn))
        end
        if iColumn == 1
            ylabel(sprintf("D = %d", iDerivative))
        end
        if iDerivative == D
            xlabel("Query point")
        end
    end
end
