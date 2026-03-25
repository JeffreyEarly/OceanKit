classdef TensorSplineUnitTests < matlab.unittest.TestCase

    methods (Test)
        function tensorInterpolatingSplineMatchesGridData(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance

            x = linspace(-1,1,6)';
            y = linspace(0,2,7)';
            [X,Y] = ndgrid(x,y);
            F = X.^2 .* Y.^3 + 2*X.*Y - 5;

            spline = InterpolatingSpline({x, y}, F, S=[3 3]);

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

            spline = InterpolatingSpline({x, y}, F, S=[3 3]);

            testCase.assertThat(spline.valueAtPoints(X, Y, D=[1 0]), IsEqualTo(dFdx, 'Within', AbsoluteTolerance(1e-9)))
        end

        function tensorSplineMatchesGriddedInterpolantSpline(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance

            x = linspace(-1,1,5)';
            y = linspace(0,2,6)';
            [X,Y] = ndgrid(x,y);
            F = sin(pi*X) + cos(0.5*pi*Y);

            spline = InterpolatingSpline({x, y}, F, S=[3 3]);
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

            spline = InterpolatingSpline({x, y}, F, S=[3 3]);
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

            spline = InterpolatingSpline({x, y}, F, S=[3 3]);

            testCase.verifyError(@() spline([X(:), Y(:)]), 'TensorSpline:InvalidEvaluationInput')
        end

        function tensorSplineValueAtPointsTreatsMatchingColumnVectorsAsPointwiseQueries(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance

            x = linspace(-1,1,6)';
            y = linspace(0,2,7)';
            [X,Y] = ndgrid(x,y);
            F = X.^2 .* Y.^3 + 2*X.*Y - 5;

            spline = InterpolatingSpline({x, y}, F, S=[3 3]);
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

            spline = InterpolatingSpline({x, y}, F, S=[3 3]);
            xq = linspace(-1, 1, 11)';
            yq = linspace(0, 2, 11)';
            expected = 2*xq.*yq.^3 + 2*yq;

            testCase.assertThat(spline.valueAtPoints(xq, yq, D=[1 0]), IsEqualTo(expected, 'Within', AbsoluteTolerance(1e-9)))
        end

        function tensorSplineValueAtPointsPreservesOneDimensionalInputShape(testCase)
            x = linspace(0, 1, 6)';
            f = (x + 1).^2;
            query = linspace(0, 1, 9);
            spline = InterpolatingSpline(x, f, S=3);

            values = spline.valueAtPoints(query);

            testCase.verifySize(values, size(query))
            testCase.verifyEqual(values, spline(query))
        end

        function tensorSplineRejectsMismatchedQueryArraySizes(testCase)
            x = linspace(-1,1,6)';
            y = linspace(0,2,7)';
            [X,Y] = ndgrid(x,y);
            F = X.^2 + Y;
            spline = InterpolatingSpline({x, y}, F, S=[3 3]);

            testCase.verifyError(@() spline.valueAtPoints((0:4)', (0:5)'), 'TensorSpline:InvalidQueryArrays')
        end

        function tensorSplineRejectsDerivativeInFunctionCallSyntax(testCase)
            x = linspace(-1,1,6)';
            y = linspace(0,2,7)';
            [X,Y] = ndgrid(x,y);
            F = X.^2 + Y;
            spline = InterpolatingSpline({x, y}, F, S=[3 3]);

            testCase.verifyError(@() spline(X, Y, [1 0]), 'TensorSpline:InvalidEvaluationInput')
        end

        function tensorSplinePlusAddsScalarOffset(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance

            x = linspace(-1,1,6)';
            y = linspace(0,2,7)';
            [X,Y] = ndgrid(x,y);
            F = X.^2 + 3*Y - 1;
            spline = InterpolatingSpline({x, y}, F, S=[3 3]);

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
            spline = InterpolatingSpline({x, y}, F, S=[3 3]);

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

            spline = InterpolatingSpline({x, y}, F, S=[3 3]);
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

            spline = InterpolatingSpline({x, y}, F, S=[1 1]);
            intspline = cumsum(spline, 2);

            testCase.assertThat(intspline(X, Y),  IsEqualTo(integralAlongY - integralAlongY(:,1), 'Within', AbsoluteTolerance(1e-10)))
        end

        function tensorSplineExposesDegreeVector(testCase)
            spline2D = InterpolatingSpline({linspace(0,1,5)', linspace(-1,1,6)'}, randn(5,6), S=[2 3]);

            testCase.verifyEqual(spline2D.S, [2 3])
        end

        function tensorSplineFevalMatchesDirectEvaluation(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance

            x = linspace(-1,1,6)';
            y = linspace(0,2,7)';
            [X,Y] = ndgrid(x,y);
            F = X.^2 .* Y.^3 + 2*X.*Y - 5;

            spline = InterpolatingSpline({x, y}, F, S=[3 3]);

            testCase.assertThat(feval(spline, X, Y), IsEqualTo(spline(X, Y), 'Within', AbsoluteTolerance(1e-10)))
            testCase.assertThat(feval(spline, X, Y, D=[1 0]), IsEqualTo(spline.valueAtPoints(X, Y, D=[1 0]), 'Within', AbsoluteTolerance(1e-10)))
        end

        function tensorSplineDomainReturnsNumericLimits(testCase)
            spline1D = InterpolatingSpline((0:4)', (0:4)', S=1);
            spline2D = InterpolatingSpline({(0:4)', (-2:2)'}, randn(5,5), S=[1 2]);

            testCase.verifyEqual(spline1D.domain, [0 4])
            testCase.verifyEqual(spline2D.domain, [0 4; -2 2])
        end

        function tensorSplineRootsRejectHigherDimensions(testCase)
            x = linspace(0,1,5)';
            y = linspace(0,1,6)';
            [X,Y] = ndgrid(x,y);
            F = X - Y;
            spline = InterpolatingSpline({x, y}, F, S=[1 1]);

            testCase.verifyError(@() roots(spline), 'TensorSpline:roots:UnsupportedDimension')
        end

        function tensorSplinePointsOfSupportMatchBasisSize(testCase)
            x = linspace(-1,1,6)';
            y = linspace(0,2,7)';
            spline = InterpolatingSpline({x, y}, x .* y', S=[3 3]);

            [supportPoints, supportVectors] = TensorSpline.pointsOfSupportFromKnotPoints(spline.knotPoints, S=spline.S);

            testCase.verifySize(supportPoints, [prod(spline.basisSize), spline.numDimensions])
            testCase.verifyEqual(cellfun(@numel, supportVectors), spline.basisSize)
        end

        function tensorSplineCoefficientsCanBeReassigned(testCase)
            tKnot = {
                [0; 0; 1; 1]
                [0; 0; 1; 1]
            };
            spline = TensorSpline(S=[1 1], knotPoints=tKnot, xi=1:4);

            spline.xi = 11:14;

            testCase.verifyEqual(spline.xi, reshape(11:14, [2 2]))
        end

        function tensorSplineRejectsInvalidCoefficientCountAssignment(testCase)
            tKnot = {
                [0; 0; 1; 1]
                [0; 0; 1; 1]
            };
            spline = TensorSpline(S=[1 1], knotPoints=tKnot, xi=1:4);

            caught = [];
            try
                spline.xi = 1:3;
            catch exception
                caught = exception;
            end

            testCase.verifyNotEmpty(caught)
            testCase.verifyEqual(caught.identifier, 'TensorSpline:InvalidCoefficientCount')
        end

        function tensorSplineOutputAffineTermsAreReadOnly(testCase)
            spline = TensorSpline(S=1, knotPoints={[0; 0; 1; 1]}, xi=[1; 2], xMean=3, xStd=4);

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
            spline = TensorSpline(S=1, knotPoints={[0; 0; 1; 1]}, xi=[1; 2]);

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
            spline = InterpolatingSpline({x, y}, F, S=[3 3]);
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
            spline = InterpolatingSpline({x, y}, F, S=[3 3]);
            rootedSpline = sqrt(spline);
            [supportPoints, supportVectors] = TensorSpline.pointsOfSupportFromKnotPoints(spline.knotPoints, S=spline.S);
            [Xsup, Ysup] = ndgrid(supportVectors{:});
            expected = reshape(TensorSpline.matrixForPointMatrix(supportPoints, knotPoints=spline.knotPoints, S=spline.S) * spline.xi(:), cellfun(@numel, supportVectors));
            expected = spline.xStd * expected + spline.xMean;

            testCase.assertThat(rootedSpline(Xsup, Ysup),  IsEqualTo(sqrt(expected), 'Within', AbsoluteTolerance(1e-10)))
        end
    end
end
