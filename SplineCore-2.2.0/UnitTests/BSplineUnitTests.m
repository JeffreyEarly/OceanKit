classdef BSplineUnitTests < matlab.unittest.TestCase

    methods (Test)
        function matrixAndPPFormsAgree(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance

            f = @(x) cos(2*pi*x/3);
            S = 3;
            t = linspace(0,3,13)';
            knotPoints = BSpline.knotPointsForDataPoints(t, S=S);

            X = BSpline.matrixForDataPoints(t, knotPoints=knotPoints, S=S);
            xi = X\f(t);
            [C,tpp] = BSpline.ppCoefficientsFromSplineCoefficients(xi=xi, knotPoints=knotPoints, S=S);

            tq = linspace(t(1),t(end),301)';
            valuesFromMatrix = BSpline.matrixForDataPoints(tq, knotPoints=knotPoints, S=S) * xi;
            valuesFromPP = BSpline.evaluateFromPPCoefficients(queryPoints=tq, C=C, tpp=tpp);

            testCase.assertThat(valuesFromPP, IsEqualTo(valuesFromMatrix, 'Within', AbsoluteTolerance(1e-10)))
        end

        function splineDOFMatchesDocumentedDataDOFMapping(testCase)
            t = linspace(0,3,13)';
            S = 3;
            oldDataDOF = 3;
            splineDOF = max(S + 1, ceil(numel(t)/oldDataDOF));
            expected = legacyKnotPointsForDataPoints(t, S + 1, oldDataDOF);
            actual = BSpline.knotPointsForDataPoints(t, S=S, splineDOF=splineDOF);

            testCase.verifyEqual(actual, expected)
        end

        function ppEvaluationHandlesUnsortedInputs(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance

            f = @(x) sin(2*pi*x/3);
            S = 3;
            t = linspace(0,3,13)';
            knotPoints = BSpline.knotPointsForDataPoints(t, S=S);

            X = BSpline.matrixForDataPoints(t, knotPoints=knotPoints, S=S);
            xi = X\f(t);
            [C,tpp] = BSpline.ppCoefficientsFromSplineCoefficients(xi=xi, knotPoints=knotPoints, S=S);

            tqSorted = linspace(t(1),t(end),31)';
            tqUnsorted = tqSorted([7 1 19 4 31 12 2 25 16 9 22 5 29 14 3 18 27 11 6 24 15 8 30 13 10 21 17 20 23 26 28]);

            expected = BSpline.matrixForDataPoints(tqUnsorted, knotPoints=knotPoints, S=S) * xi;
            actual = BSpline.evaluateFromPPCoefficients(queryPoints=tqUnsorted, C=C, tpp=tpp);

            testCase.assertThat(actual, IsEqualTo(expected, 'Within', AbsoluteTolerance(1e-10)))
        end

        function fevalDelegatesToValueAtPoints(testCase)
            t = linspace(0,1,9)';
            x = sin(2*pi*t);
            knotPoints = BSpline.knotPointsForDataPoints(t, S=3);
            spline = BSpline(S=3, knotPoints=knotPoints, xi=BSpline.matrixForDataPoints(t, knotPoints=knotPoints, S=3)\x);
            tq = linspace(0,1,21)';

            testCase.verifyEqual(feval(spline, tq), spline.valueAtPoints(tq))
            testCase.verifyEqual(feval(spline, tq, D=1), spline.valueAtPoints(tq, D=1))
        end

        function clearingCoefficientsClearsCachedState(testCase)
            K = 3;
            S = K - 1;
            t = linspace(0,1,5)';
            knotPoints = [t(1)*ones(K,1); t(2:end-1); t(end)*ones(K,1)];
            spline = BSpline(S=S, knotPoints=knotPoints, xi=ones(length(knotPoints)-K,1));

            spline.xi = [];

            testCase.verifyEmpty(spline.C)
            testCase.verifyEmpty(spline.t_pp)
            testCase.verifyEmpty(spline.Xtpp)
        end

        function settingCoefficientsRebuildsCachedState(testCase)
            K = 3;
            S = K - 1;
            t = linspace(0,1,5)';
            knotPoints = [t(1)*ones(K,1); t(2:end-1); t(end)*ones(K,1)];
            spline = BSpline(S=S, knotPoints=knotPoints, xi=[]);

            spline.xi = ones(length(knotPoints)-K,1);

            testCase.verifyNotEmpty(spline.C)
            testCase.verifyNotEmpty(spline.t_pp)
            testCase.verifyNotEmpty(spline.Xtpp)
        end

        function knotSequenceIsReadOnly(testCase)
            K = 3;
            S = K - 1;
            t = linspace(0,1,5)';
            knotPoints = [t(1)*ones(K,1); t(2:end-1); t(end)*ones(K,1)];
            spline = BSpline(S=S, knotPoints=knotPoints, xi=ones(length(knotPoints)-K,1));

            caught = [];
            try
                spline.knotPoints = [0; 0; 0; 0.5; 1; 1; 1];
            catch exception
                caught = exception;
            end

            testCase.verifyNotEmpty(caught)
            testCase.verifyEqual(spline.knotPoints, knotPoints)
            testCase.verifyEqual(spline.xi, ones(length(knotPoints)-K,1))
        end

        function orderIsReadOnly(testCase)
            spline = BSpline(S=2, knotPoints=[0; 0; 0; 1; 1; 1], xi=[1; 2; 3]);

            caught = [];
            try
                spline.K = 4;
            catch exception
                caught = exception;
            end

            testCase.verifyNotEmpty(caught)
            testCase.verifyEqual(spline.S, 2)
            testCase.verifyEqual(spline.K, 3)
        end

        function outputAffineTermsAreReadOnly(testCase)
            spline = BSpline(S=2, knotPoints=[0; 0; 0; 1; 1; 1], xi=[1; 2; 3], xMean=4, xStd=5);

            xMeanException = [];
            try
                spline.xMean = 7;
            catch exception
                xMeanException = exception;
            end

            xStdException = [];
            try
                spline.xStd = 8;
            catch exception
                xStdException = exception;
            end

            testCase.verifyNotEmpty(xMeanException)
            testCase.verifyNotEmpty(xStdException)
            testCase.verifyEqual(spline.xMean, 4)
            testCase.verifyEqual(spline.xStd, 5)
        end

        function integratedSplineStateMatchesLegacyCumsumWithAffineTerms(testCase)
            t = linspace(-1, 1, 9)';
            S = 3;
            knotPoints = BSpline.knotPointsForDataPoints(t, S=S);
            xi = reshape(linspace(-0.6, 0.9, numel(knotPoints) - S - 1), [], 1);
            xMean = 2.5;
            xStd = 1.75;

            [xiIntegrated, knotPointsIntegrated, SIntegrated] = BSpline.integratedSplineState( ...
                xi, knotPoints=knotPoints, S=S, xMean=xMean, xStd=xStd);
            [xiExpected, knotPointsExpected, SExpected] = legacyIntegratedSplineState(xi, knotPoints, S, xMean, xStd);

            testCase.verifyEqual(xiIntegrated, xiExpected, AbsTol=1e-12)
            testCase.verifyEqual(knotPointsIntegrated, knotPointsExpected, AbsTol=1e-12)
            testCase.verifyEqual(SIntegrated, SExpected)

            referenceSpline = BSpline(S=S, knotPoints=knotPoints, xi=xi, xMean=xMean, xStd=xStd);
            integratedSpline = BSpline(S=SIntegrated, knotPoints=knotPointsIntegrated, xi=xiIntegrated);
            queryPoints = linspace(t(1), t(end), 41)';
            referenceIntegratedSpline = cumsum(referenceSpline);

            testCase.verifyEqual(integratedSpline(queryPoints), referenceIntegratedSpline(queryPoints), AbsTol=1e-12)
        end

        function integratedSplineStateMatchesLegacyColumnwiseForMatrixCoefficients(testCase)
            t = linspace(0, 1, 7)';
            S = 2;
            knotPoints = BSpline.knotPointsForDataPoints(t, S=S);
            numCoefficients = numel(knotPoints) - S - 1;
            xi = reshape(linspace(-1.2, 1.1, 3 * numCoefficients), numCoefficients, 3);
            xMean = -0.75;
            xStd = 0.5;

            [xiIntegrated, knotPointsIntegrated, SIntegrated] = BSpline.integratedSplineState( ...
                xi, knotPoints=knotPoints, S=S, xMean=xMean, xStd=xStd);
            [xiExpected, knotPointsExpected, SExpected] = legacyIntegratedSplineState(xi, knotPoints, S, xMean, xStd);

            testCase.verifyEqual(xiIntegrated, xiExpected, AbsTol=1e-12)
            testCase.verifyEqual(knotPointsIntegrated, knotPointsExpected, AbsTol=1e-12)
            testCase.verifyEqual(SIntegrated, SExpected)
        end

        function integralMatrixForDataPointsMatchesInterpolatingSplineAntiderivative(testCase)
            t = sort([0; 0.1; 0.25; 0.4; 0.7; 0.9; 1.0]);
            values = cos(2 * pi * t) + 0.2 * t;
            knotPoints = BSpline.knotPointsForDataPoints(t, S=3);
            spline = InterpolatingSpline.fromGriddedValues(t, values, S=3);
            integratedSpline = cumsum(spline);

            scalarWeights = BSpline.integralMatrixForDataPoints(t, t(end), knotPoints=knotPoints, S=3);
            vectorQueryPoints = linspace(t(1), t(end), 17)';
            vectorWeights = BSpline.integralMatrixForDataPoints(t, vectorQueryPoints, knotPoints=knotPoints, S=3);

            testCase.verifyEqual(scalarWeights * values, integratedSpline(t(end)), AbsTol=1e-12)
            testCase.verifyEqual(vectorWeights * values, integratedSpline(vectorQueryPoints), AbsTol=1e-12)
        end

        function powerPreservesPositivityForNearlyNonnegativeSupportValues(testCase)
            t = linspace(0,1,9)';
            x = [0.3214; 0.0927; 0.1944; 0.0228; 0.0278; 0.5334; 0.9278; 0.0508; 0.3518];
            exponent = 1.376389;
            knotPoints = BSpline.knotPointsForDataPoints(t, S=3);
            spline = BSpline(S=3, knotPoints=knotPoints, xi=BSpline.matrixForDataPoints(t, knotPoints=knotPoints, S=3)\x);

            poweredSpline = spline.^exponent;

            supportPoints = BSpline.pointsOfSupportFromKnotPoints(spline.knotPoints, S=spline.S);
            supportValues = spline.valueAtPoints(supportPoints);
            supportValues(abs(supportValues) < 2*eps) = 0;
            poweredK = ceil(exponent*spline.K);
            poweredKnotPoints = BSpline.knotPointsForDataPoints(supportPoints, S=poweredK-1);
            X = BSpline.matrixForDataPoints(supportPoints, knotPoints=poweredKnotPoints, S=poweredK-1);
            unconstrainedSpline = BSpline(S=poweredK-1, knotPoints=poweredKnotPoints, xi=X\(supportValues.^exponent));
            tq = linspace(0,1,1001)';

            testCase.verifyLessThan(min(unconstrainedSpline(tq)), -1e-3)
            testCase.verifyGreaterThanOrEqual(min(poweredSpline(tq)), -1e-10)
        end
    end

end

function [xiIntegrated, knotPointsIntegrated, SIntegrated] = legacyIntegratedSplineState(xi, knotPoints, S, xMean, xStd)
xi = reshape(xi, size(xi, 1), []);
if abs(xMean) > 0 || abs(xStd - 1) > 0
    supportPoints = BSpline.pointsOfSupportFromKnotPoints(knotPoints, S=S);
    basisMatrix = BSpline.matrixForDataPoints(supportPoints, knotPoints=knotPoints, S=S);
    xi = xStd * xi + basisMatrix \ (xMean * ones(numel(supportPoints), 1));
end

K = S + 1;
numCoefficients = size(xi, 1);
dt = (knotPoints(1 + K:numCoefficients + K) - knotPoints(1:numCoefficients)) / K;
xiIntegrated = [zeros(1, size(xi, 2), 'like', xi); cumsum(xi .* reshape(dt, [], 1), 1)];
knotPointsIntegrated = [knotPoints(1); knotPoints; knotPoints(end)];
SIntegrated = S + 1;
end

function tKnot = legacyKnotPointsForDataPoints(t, K, dataDOF)
tData = sort(t);
tData = [tData(1); tData(1+dataDOF:dataDOF:end-dataDOF); tData(end)];
mustBeGreaterThanOrEqual(numel(tData), K);

tPseudo = interp1((0:numel(tData)-1)', tData, linspace(0, numel(tData)-1, numel(tData)).');
if mod(K, 2) == 1
    dt = diff(tPseudo);
    tKnot = [tPseudo(1); tPseudo(1:end-1) + dt/2; tPseudo(end)];
    for i = 1:((K-1)/2)
        tKnot(2) = [];
        tKnot(end-1) = [];
    end
else
    tKnot = tPseudo;
    for i = 1:((K-2)/2)
        tKnot(2) = [];
        tKnot(end-1) = [];
    end
end

tKnot = [repmat(tKnot(1), K-1, 1); tKnot; repmat(tKnot(end), K-1, 1)];
end
