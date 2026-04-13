classdef TensorSplineUnitTests < matlab.unittest.TestCase

    methods (Test)
        function tensorInterpolatingSplineMatchesGridData(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance

            x = linspace(-1,1,6)';
            y = linspace(0,2,7)';
            [X,Y] = ndgrid(x,y);
            F = X.^2 .* Y.^3 + 2*X.*Y - 5;

            spline = InterpolatingSpline.fromGriddedValues({x, y}, F, S=[3 3]);

            testCase.assertThat(spline(X, Y), IsEqualTo(F, 'Within', AbsoluteTolerance(1e-10)))
        end

        function tensorSplineMatchesPolynomialDerivative(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance

            x = linspace(-1,1,6)';
            y = linspace(0,2,7)';
            [X,Y] = ndgrid(x,y);
            F = X.^2 .* Y.^3 + 2*X.*Y - 5;
            dFdx = 2*X.*Y.^3 + 2*Y;

            spline = InterpolatingSpline.fromGriddedValues({x, y}, F, S=[3 3]);

            testCase.assertThat(spline.valueAtPoints(X, Y, D=[1 0]), IsEqualTo(dFdx, 'Within', AbsoluteTolerance(1e-9)))
        end

        function tensorSplineMatchesGriddedInterpolantSpline(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance

            x = linspace(-1,1,5)';
            y = linspace(0,2,6)';
            [X,Y] = ndgrid(x,y);
            F = sin(pi*X) + cos(0.5*pi*Y);

            spline = InterpolatingSpline.fromGriddedValues({x, y}, F, S=[3 3]);
            interpolant = griddedInterpolant({x,y}, F, 'spline');

            xq = linspace(x(1),x(end),21)';
            yq = linspace(y(1),y(end),25)';
            [Xq,Yq] = ndgrid(xq,yq);

            expected = interpolant({xq,yq});
            actual = spline(Xq, Yq);

            testCase.assertThat(actual, IsEqualTo(expected, 'Within', AbsoluteTolerance(1e-10)))
        end

        function tensorSplineMatrixMatchesEvaluation(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance

            x = linspace(-1,1,6)';
            y = linspace(0,2,7)';
            [X,Y] = ndgrid(x,y);
            F = X.^2 .* Y.^3 + 2*X.*Y - 5;

            spline = InterpolatingSpline.fromGriddedValues({x, y}, F, S=[3 3]);
            queryPoints = [X(:), Y(:)];
            basisMatrix = TensorSpline.matrixForPointMatrix(queryPoints, knotPoints=spline.knotPoints, S=spline.S);
            values = reshape(basisMatrix * spline.xi(:), size(F));
            values = spline.xStd * values + spline.xMean;

            testCase.assertThat(values, IsEqualTo(F, 'Within', AbsoluteTolerance(1e-10)))
        end

        function tensorSplineRejectsPointMatrixEvaluationSyntax(testCase)
            x = linspace(-1,1,6)';
            y = linspace(0,2,7)';
            [X,Y] = ndgrid(x,y);
            F = X.^2 + Y;

            spline = InterpolatingSpline.fromGriddedValues({x, y}, F, S=[3 3]);

            testCase.verifyError(@() spline([X(:), Y(:)]), 'TensorSpline:InvalidEvaluationInput')
        end

        function tensorSplineValueAtPointsTreatsMatchingColumnVectorsAsPointwiseQueries(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance

            x = linspace(-1,1,6)';
            y = linspace(0,2,7)';
            [X,Y] = ndgrid(x,y);
            F = X.^2 .* Y.^3 + 2*X.*Y - 5;

            spline = InterpolatingSpline.fromGriddedValues({x, y}, F, S=[3 3]);
            xq = linspace(-1, 1, 11)';
            yq = linspace(0, 2, 11)';
            expected = xq.^2 .* yq.^3 + 2*xq.*yq - 5;

            testCase.assertThat(spline.valueAtPoints(xq, yq), IsEqualTo(expected, 'Within', AbsoluteTolerance(1e-10)))
            testCase.assertThat(spline(xq, yq), IsEqualTo(expected, 'Within', AbsoluteTolerance(1e-10)))
        end

        function tensorSplineValueAtPointsSupportsDerivativeOrders(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance

            x = linspace(-1,1,6)';
            y = linspace(0,2,7)';
            [X,Y] = ndgrid(x,y);
            F = X.^2 .* Y.^3 + 2*X.*Y - 5;
            dFdx = 2*X.*Y.^3 + 2*Y;

            spline = InterpolatingSpline.fromGriddedValues({x, y}, F, S=[3 3]);
            xq = linspace(-1, 1, 11)';
            yq = linspace(0, 2, 11)';
            expected = 2*xq.*yq.^3 + 2*yq;

            testCase.assertThat(spline.valueAtPoints(xq, yq, D=[1 0]), IsEqualTo(expected, 'Within', AbsoluteTolerance(1e-9)))
        end

        function tensorSplineValueAtPointsPreservesOneDimensionalInputShape(testCase)
            x = linspace(0, 1, 6)';
            f = (x + 1).^2;
            query = linspace(0, 1, 9);
            spline = InterpolatingSpline.fromGriddedValues(x, f, S=3);

            values = spline.valueAtPoints(query);

            testCase.verifySize(values, size(query))
            testCase.verifyEqual(values, spline(query))
        end

        function tensorSplineRejectsMismatchedQueryArraySizes(testCase)
            x = linspace(-1,1,6)';
            y = linspace(0,2,7)';
            [X,Y] = ndgrid(x,y);
            F = X.^2 + Y;
            spline = InterpolatingSpline.fromGriddedValues({x, y}, F, S=[3 3]);

            testCase.verifyError(@() spline.valueAtPoints((0:4)', (0:5)'), 'TensorSpline:InvalidQueryArrays')
        end

        function tensorSplineRejectsDerivativeInFunctionCallSyntax(testCase)
            x = linspace(-1,1,6)';
            y = linspace(0,2,7)';
            [X,Y] = ndgrid(x,y);
            F = X.^2 + Y;
            spline = InterpolatingSpline.fromGriddedValues({x, y}, F, S=[3 3]);

            testCase.verifyError(@() spline(X, Y, [1 0]), 'TensorSpline:InvalidEvaluationInput')
        end

        function tensorSplinePlusAddsScalarOffset(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance

            x = linspace(-1,1,6)';
            y = linspace(0,2,7)';
            [X,Y] = ndgrid(x,y);
            F = X.^2 + 3*Y - 1;
            spline = InterpolatingSpline.fromGriddedValues({x, y}, F, S=[3 3]);

            shifted = spline + 1;

            testCase.assertThat(shifted(X, Y), IsEqualTo(F + 1, 'Within', AbsoluteTolerance(1e-10)))
        end

        function tensorSplineMtimesScalesSpline(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance

            x = linspace(-1,1,6)';
            y = linspace(0,2,7)';
            [X,Y] = ndgrid(x,y);
            F = X.^2 + 3*Y - 1;
            spline = InterpolatingSpline.fromGriddedValues({x, y}, F, S=[3 3]);

            scaled = -2 * spline;

            testCase.assertThat(scaled(X, Y), IsEqualTo(-2*F, 'Within', AbsoluteTolerance(1e-10)))
        end

        function tensorSplineDiffReturnsMixedPartialSpline(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance

            x = linspace(-1,1,6)';
            y = linspace(0,2,7)';
            [X,Y] = ndgrid(x,y);
            F = X.^2 .* Y.^3 + 2*X.*Y - 5;
            d2Fdxy = 6*X.*Y.^2 + 2;

            spline = InterpolatingSpline.fromGriddedValues({x, y}, F, S=[3 3]);
            dspline = diff(spline, [1 1]);

            testCase.assertThat(dspline(X, Y), IsEqualTo(d2Fdxy, 'Within', AbsoluteTolerance(1e-9)))
        end

        function tensorSplineCumsumMatchesIntegralAlongSpecifiedDimension(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance

            x = linspace(-1,1,6)';
            y = linspace(0,2,7)';
            [X,Y] = ndgrid(x,y);
            F = X + 2*Y;
            integralAlongY = X.*Y + Y.^2;

            spline = InterpolatingSpline.fromGriddedValues({x, y}, F, S=[1 1]);
            intspline = cumsum(spline, 2);

            testCase.assertThat(intspline(X, Y),  IsEqualTo(integralAlongY - integralAlongY(:,1), 'Within', AbsoluteTolerance(1e-10)))
        end

        function tensorSplineCumsumWithAffineTermsMatchesLegacySliceIntegration(testCase)
            x = linspace(-1, 1, 6)';
            y = linspace(0, 2, 7)';
            knotPoints = {BSpline.knotPointsForDataPoints(x, S=1), BSpline.knotPointsForDataPoints(y, S=1)};
            xi = reshape(linspace(-0.4, 0.6, 42), 6, 7);
            spline = TensorSpline(S=[1 1], knotAxes=SplineAxis.arrayFromVectors(knotPoints), xi=xi, xMean=3, xStd=4);
            splineKnotPoints = spline.knotPoints;

            intspline = cumsum(spline, 2);
            [expectedXi, expectedTKnotDim, expectedK] = legacyIntegrateAlongDimension( ...
                spline.xi, splineKnotPoints{2}, spline.K(2), 2, spline.xMean, spline.xStd);
            expectedSpline = TensorSpline(S=[spline.S(1), expectedK - 1], ...
                knotAxes=SplineAxis.arrayFromVectors({splineKnotPoints{1}, expectedTKnotDim}), xi=expectedXi);
            [X, Y] = ndgrid(linspace(-1, 1, 11)', linspace(0, 2, 13)');

            testCase.verifyEqual(intspline(X, Y), expectedSpline(X, Y), AbsTol=1e-12)
        end

        function tensorSplineExposesDegreeVector(testCase)
            spline2D = InterpolatingSpline.fromGriddedValues({linspace(0,1,5)', linspace(-1,1,6)'}, randn(5,6), S=[2 3]);

            testCase.verifyEqual(spline2D.S, [2 3])
        end

        function tensorSplineFevalMatchesDirectEvaluation(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance

            x = linspace(-1,1,6)';
            y = linspace(0,2,7)';
            [X,Y] = ndgrid(x,y);
            F = X.^2 .* Y.^3 + 2*X.*Y - 5;

            spline = InterpolatingSpline.fromGriddedValues({x, y}, F, S=[3 3]);

            testCase.assertThat(feval(spline, X, Y), IsEqualTo(spline(X, Y), 'Within', AbsoluteTolerance(1e-10)))
            testCase.assertThat(feval(spline, X, Y, D=[1 0]), IsEqualTo(spline.valueAtPoints(X, Y, D=[1 0]), 'Within', AbsoluteTolerance(1e-10)))
        end

        function tensorSplineDomainReturnsNumericLimits(testCase)
            spline1D = InterpolatingSpline.fromGriddedValues((0:4)', (0:4)', S=1);
            spline2D = InterpolatingSpline.fromGriddedValues({(0:4)', (-2:2)'}, randn(5,5), S=[1 2]);

            testCase.verifyEqual(spline1D.domain, [0 4])
            testCase.verifyEqual(spline2D.domain, [0 4; -2 2])
        end

        function tensorSplineRootsRejectHigherDimensions(testCase)
            x = linspace(0,1,5)';
            y = linspace(0,1,6)';
            [X,Y] = ndgrid(x,y);
            F = X - Y;
            spline = InterpolatingSpline.fromGriddedValues({x, y}, F, S=[1 1]);

            testCase.verifyError(@() roots(spline), 'TensorSpline:roots:UnsupportedDimension')
        end

        function tensorSplinePointsOfSupportMatchBasisSize(testCase)
            x = linspace(-1,1,6)';
            y = linspace(0,2,7)';
            spline = InterpolatingSpline.fromGriddedValues({x, y}, x .* y', S=[3 3]);

            [supportPoints, supportVectors] = TensorSpline.pointsOfSupportFromKnotPoints(spline.knotPoints, S=spline.S);

            testCase.verifySize(supportPoints, [prod(spline.basisSize), spline.numDimensions])
            testCase.verifyEqual(cellfun(@numel, supportVectors), spline.basisSize)
        end

        function tensorSplineCoefficientsCanBeReassigned(testCase)
            tKnot = {
                [0; 0; 1; 1]
                [0; 0; 1; 1]
            };
            spline = TensorSpline.fromKnotPoints(tKnot, 1:4, S=[1 1]);

            spline.xi = 11:14;

            testCase.verifyEqual(spline.xi, reshape(11:14, [2 2]))
        end

        function tensorSplineRejectsInvalidCoefficientCountAssignment(testCase)
            tKnot = {
                [0; 0; 1; 1]
                [0; 0; 1; 1]
            };
            spline = TensorSpline.fromKnotPoints(tKnot, 1:4, S=[1 1]);

            caught = [];
            try
                spline.xi = 1:3;
            catch exception
                caught = exception;
            end

            testCase.verifyNotEmpty(caught)
            testCase.verifyEqual(caught.identifier, 'TensorSpline:InvalidCoefficientCount')
        end

        function tensorSplineLowLevelConstructorAcceptsPublicAxisAPI(testCase)
            knotPoints = {
                [0; 0; 1; 1]
                [-1; -1; 1; 1]
            };
            xi = 1:4;

            spline = TensorSpline(S=[1 1], knotAxes=SplineAxis.arrayFromVectors(knotPoints), xi=xi, xMean=3, xStd=4);
            knotAxes = spline.knotAxes;
            returnedKnotPoints = spline.knotPoints;

            testCase.verifyEqual(numel(knotAxes), 2)
            testCase.verifyEqual(knotAxes(1).values, knotPoints{1})
            testCase.verifyEqual(knotAxes(2).values, knotPoints{2})
            testCase.verifyEqual(returnedKnotPoints{1}, knotPoints{1})
            testCase.verifyEqual(returnedKnotPoints{2}, knotPoints{2})
            testCase.verifyEqual(spline.xi, reshape(xi, [2 2]))
            testCase.verifyEqual(spline.xMean, 3)
            testCase.verifyEqual(spline.xStd, 4)
        end

        function tensorSplineLowLevelConstructorRejectsMissingKnotAxes(testCase)
            caught = [];
            try
                TensorSpline(S=1, knotAxes=SplineAxis.empty(0,1), xi=[1; 2]);
            catch exception
                caught = exception;
            end

            testCase.verifyNotEmpty(caught)
        end

        function tensorSplineOutputAffineTermsAreReadOnly(testCase)
            spline = TensorSpline.fromKnotPoints([0; 0; 1; 1], [1; 2], S=1, xMean=3, xStd=4);

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
            testCase.verifyEqual(spline.xMean, 3)
            testCase.verifyEqual(spline.xStd, 4)
        end

        function tensorSplineOrderIsReadOnly(testCase)
            spline = TensorSpline.fromKnotPoints([0; 0; 1; 1], [1; 2], S=1);

            caught = [];
            try
                spline.K = 4;
            catch exception
                caught = exception;
            end

            testCase.verifyNotEmpty(caught)
            testCase.verifyEqual(spline.S, 1)
            testCase.verifyEqual(spline.K, 2)
        end

        function tensorSplinePowerMatchesSquaredValuesOnSupport(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance

            x = linspace(-1,1,6)';
            y = linspace(0,2,7)';
            [X,Y] = ndgrid(x,y);
            F = X.^2 + Y + 2;
            spline = InterpolatingSpline.fromGriddedValues({x, y}, F, S=[3 3]);
            squaredSpline = spline.^2;
            [supportPoints, supportVectors] = TensorSpline.pointsOfSupportFromKnotPoints(spline.knotPoints, S=spline.S);
            [Xsup, Ysup] = ndgrid(supportVectors{:});
            expected = reshape(TensorSpline.matrixForPointMatrix(supportPoints, knotPoints=spline.knotPoints, S=spline.S) * spline.xi(:), cellfun(@numel, supportVectors));
            expected = spline.xStd * expected + spline.xMean;

            testCase.assertThat(squaredSpline(Xsup, Ysup),  IsEqualTo(expected.^2, 'Within', AbsoluteTolerance(1e-10)))
        end

        function tensorSplineSqrtMatchesSquareRootValuesOnSupport(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance

            x = linspace(0,1,6)';
            y = linspace(0,2,7)';
            [X,Y] = ndgrid(x,y);
            F = (X + 1).^2 + (Y + 2).^2;
            spline = InterpolatingSpline.fromGriddedValues({x, y}, F, S=[3 3]);
            rootedSpline = sqrt(spline);
            [supportPoints, supportVectors] = TensorSpline.pointsOfSupportFromKnotPoints(spline.knotPoints, S=spline.S);
            [Xsup, Ysup] = ndgrid(supportVectors{:});
            expected = reshape(TensorSpline.matrixForPointMatrix(supportPoints, knotPoints=spline.knotPoints, S=spline.S) * spline.xi(:), cellfun(@numel, supportVectors));
            expected = spline.xStd * expected + spline.xMean;

            testCase.assertThat(rootedSpline(Xsup, Ysup),  IsEqualTo(sqrt(expected), 'Within', AbsoluteTolerance(1e-10)))
        end
    end
end

function [xi, tKnot, K] = legacyIntegrateAlongDimension(xi, tKnot, K, dim, xMean, xStd)
perm = [dim, 1:(dim-1), (dim+1):ndims(xi)];
xiPermuted = permute(xi, perm);
xiMatrix = reshape(xiPermuted, size(xiPermuted, 1), []);
[transformedMatrix, tKnot, SIntegrated] = legacyIntegratedSplineState(xiMatrix, tKnot, K - 1, xMean, xStd);
K = SIntegrated + 1;

outputSize = size(xiPermuted);
outputSize(1) = size(transformedMatrix, 1);
xi = ipermute(reshape(transformedMatrix, outputSize), perm);
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
