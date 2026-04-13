classdef TrajectorySplineUnitTests < matlab.unittest.TestCase

    methods (Test)
        function trajectorySplineConstructsFromValidSamples(testCase)
            t = linspace(0, 1, 11)';
            x = cos(2*pi*t);
            y = sin(2*pi*t);

            trajectory = TrajectorySpline.fromData(t, x, y, S=3);

            testCase.verifyClass(trajectory, 'TrajectorySpline')
            testCase.verifyClass(trajectory.x, 'ConstrainedSpline')
            testCase.verifyClass(trajectory.y, 'ConstrainedSpline')
            testCase.verifyEqual(trajectory.t, t)
        end

        function trajectorySplineInterpolatesXSamples(testCase)
            import matlab.unittest.constraints.AbsoluteTolerance
            import matlab.unittest.constraints.IsEqualTo

            t = linspace(0, 1, 11)';
            x = cos(2*pi*t);
            y = sin(2*pi*t);
            trajectory = TrajectorySpline.fromData(t, x, y, S=3);

            testCase.assertThat(trajectory.x(t), IsEqualTo(x, 'Within', AbsoluteTolerance(10*eps)))
        end

        function trajectorySplineInterpolatesYSamples(testCase)
            import matlab.unittest.constraints.AbsoluteTolerance
            import matlab.unittest.constraints.IsEqualTo

            t = linspace(0, 1, 11)';
            x = cos(2*pi*t);
            y = sin(2*pi*t);
            trajectory = TrajectorySpline.fromData(t, x, y, S=3);

            testCase.assertThat(trajectory.y(t), IsEqualTo(y, 'Within', AbsoluteTolerance(10*eps)))
        end

        function trajectorySplinePropagatesSharedDegree(testCase)
            t = linspace(0, 1, 11)';
            x = cos(2*pi*t);
            y = sin(2*pi*t);

            trajectory = TrajectorySpline.fromData(t, x, y, S=5);

            testCase.verifyEqual(trajectory.x.S, 5)
            testCase.verifyEqual(trajectory.y.S, 5)
            testCase.verifyEqual(trajectory.x.K, 6)
            testCase.verifyEqual(trajectory.y.K, 6)
        end

        function trajectorySplineExposesComponentFitData(testCase)
            t = linspace(0, 1, 11)';
            x = cos(2*pi*t);
            y = sin(2*pi*t);

            trajectory = TrajectorySpline.fromData(t, x, y, S=3);

            testCase.verifyEqual(trajectory.x.dataPoints, t)
            testCase.verifyEqual(trajectory.y.dataPoints, t)
            testCase.verifyEqual(trajectory.x.dataValues, x)
            testCase.verifyEqual(trajectory.y.dataValues, y)
        end

        function trajectorySplineEvaluatesUAsXDerivative(testCase)
            t = linspace(0, 1, 11)';
            x = cos(2*pi*t);
            y = sin(2*pi*t);
            trajectory = TrajectorySpline.fromData(t, x, y, S=3);

            testCase.verifyEqual(trajectory.u(t), trajectory.x.valueAtPoints(t, D=1), AbsTol=10*eps)
        end

        function trajectorySplineEvaluatesVAsYDerivative(testCase)
            t = linspace(0, 1, 11)';
            x = cos(2*pi*t);
            y = sin(2*pi*t);
            trajectory = TrajectorySpline.fromData(t, x, y, S=3);

            testCase.verifyEqual(trajectory.v(t), trajectory.y.valueAtPoints(t, D=1), AbsTol=10*eps)
        end

        function trajectorySplineVelocityMethodsPreserveQueryShape(testCase)
            t = linspace(0, 1, 11)';
            x = cos(2*pi*t);
            y = sin(2*pi*t);
            trajectory = TrajectorySpline.fromData(t, x, y, S=3);

            rowQuery = reshape(linspace(0.1, 0.9, 5), 1, []);
            gridQuery = reshape(linspace(0.1, 0.9, 6), 2, 3);

            testCase.verifySize(trajectory.u(rowQuery), size(rowQuery))
            testCase.verifySize(trajectory.v(rowQuery), size(rowQuery))
            testCase.verifySize(trajectory.u(gridQuery), size(gridQuery))
            testCase.verifySize(trajectory.v(gridQuery), size(gridQuery))
        end

        function trajectorySplineFromDataStoresParameterAsColumnVector(testCase)
            t = linspace(0, 1, 11);
            x = cos(2*pi*t);
            y = sin(2*pi*t);

            trajectory = TrajectorySpline.fromData(t, x, y, S=3);

            testCase.verifySize(trajectory.t, [numel(t), 1])
            testCase.verifyEqual(trajectory.t, t(:))
        end

        function trajectorySplineRejectsMismatchedVectorLengths(testCase)
            t = linspace(0, 1, 11)';
            x = cos(2*pi*t);
            y = sin(2*pi*t(1:end-1));

            testCase.verifyError(@() TrajectorySpline.fromData(t, x, y, S=3), 'TrajectorySpline:SizeMismatch')
        end

        function trajectorySplineFromDataRejectsEmptyInputs(testCase)
            testCase.verifyError(@() TrajectorySpline.fromData([], [], [], S=3), 'MATLAB:validators:mustBeNonempty')
        end

        function trajectorySplineFromDataRejectsNonmonotonicParameter(testCase)
            t = [0; 0.5; 0.5; 1];
            x = [1; 0; -1; 0];
            y = [0; 1; 0; -1];

            testCase.verifyError(@() TrajectorySpline.fromData(t, x, y, S=3), 'TrajectorySpline:NonmonotonicParameter')
        end

        function trajectorySplineConstructsFromCanonicalComponentSplines(testCase)
            t = linspace(0, 1, 11);
            x = cos(2*pi*t);
            y = sin(2*pi*t);
            xSpline = ConstrainedSpline.fromData(t', x', S=3);
            ySpline = ConstrainedSpline.fromData(t', y', S=3);

            trajectory = TrajectorySpline(t=t, x=xSpline, y=ySpline);

            testCase.verifyEqual(trajectory.t, t(:))
            testCase.verifySameHandle(trajectory.x, xSpline)
            testCase.verifySameHandle(trajectory.y, ySpline)
            testCase.verifyEqual(trajectory.x(t'), x(:), AbsTol=10*eps)
            testCase.verifyEqual(trajectory.y(t'), y(:), AbsTol=10*eps)
            testCase.verifyEqual(trajectory.u(t'), xSpline.valueAtPoints(t', D=1), AbsTol=10*eps)
            testCase.verifyEqual(trajectory.v(t'), ySpline.valueAtPoints(t', D=1), AbsTol=10*eps)
        end

        function trajectorySplineCanonicalConstructorRejectsNonmonotonicParameter(testCase)
            t = [0; 0.5; 0.5; 1];
            xSpline = ConstrainedSpline.fromData([0; 0.5; 1], [0; 1; 0], S=2);
            ySpline = ConstrainedSpline.fromData([0; 0.5; 1], [0; 0; 1], S=2);

            testCase.verifyError(@() TrajectorySpline(t=t, x=xSpline, y=ySpline), 'TrajectorySpline:NonmonotonicParameter')
        end

        function trajectorySplineCanonicalConstructorRejectsMultidimensionalSplines(testCase)
            t = linspace(0, 1, 5)';
            xSpline = TensorSpline.fromKnotPoints({[-1; -1; 1; 1], [-1; -1; 1; 1]}, zeros(4, 1), S=[1 1]);
            ySpline = ConstrainedSpline.fromData(t, sin(2*pi*t), S=3);

            testCase.verifyError(@() TrajectorySpline(t=t, x=xSpline, y=ySpline), 'TrajectorySpline:InvalidComponentSpline')
        end
    end
end
