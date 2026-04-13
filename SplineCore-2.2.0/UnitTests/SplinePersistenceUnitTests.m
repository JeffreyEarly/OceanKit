classdef SplinePersistenceUnitTests < matlab.unittest.TestCase

    methods (Test)
        function bSplineRoundTripPersistsOnlyMinimalState(testCase)
            t = linspace(0, 1, 11)';
            knotPoints = BSpline.knotPointsForDataPoints(t, S=3);
            values = sin(2*pi*t);
            X = BSpline.matrixForDataPoints(t, knotPoints=knotPoints, S=3);
            spline = BSpline(S=3, knotPoints=knotPoints, xi=X\values, xMean=2, xStd=0.5);

            path = string(strcat(tempname, '.nc'));
            cleanupPath = onCleanup(@() deleteIfPresent(path)); %#ok<NASGU>
            spline.writeToFile(path, shouldOverwriteExisting=true);

            ncfile = NetCDFFile(path, shouldReadOnly=true);
            cleanupFile = onCleanup(@() ncfile.close()); %#ok<NASGU>
            loaded = BSpline.annotatedClassFromFile(path);
            tq = linspace(0, 1, 41)';

            testCase.verifyFalse(ncfile.hasVariableWithName('Xtpp'))
            testCase.verifyFalse(ncfile.hasVariableWithName('t_pp'))
            testCase.verifyFalse(ncfile.hasVariableWithName('C'))
            testCase.verifyEqual(loaded.S, spline.S)
            testCase.verifyEqual(loaded.knotPoints, spline.knotPoints)
            testCase.verifyEqual(loaded.xi, spline.xi, AbsTol=10*eps)
            testCase.verifyEqual(loaded.xMean, spline.xMean)
            testCase.verifyEqual(loaded.xStd, spline.xStd)
            testCase.verifyEqual(loaded(tq), spline(tq), AbsTol=1e-12)
        end

        function tensorSplineOneDimensionalRoundTripUsesKnotAxes(testCase)
            t = linspace(-1, 1, 9)';
            knotPoints = BSpline.knotPointsForDataPoints(t, S=3);
            xi = (1:(numel(knotPoints) - 4)).';
            spline = TensorSpline.fromKnotPoints(knotPoints, xi, S=3, xMean=1.5, xStd=0.25);

            path = string(strcat(tempname, '.nc'));
            cleanupPath = onCleanup(@() deleteIfPresent(path)); %#ok<NASGU>
            spline.writeToFile(path, shouldOverwriteExisting=true);

            ncfile = NetCDFFile(path, shouldReadOnly=true);
            cleanupFile = onCleanup(@() ncfile.close()); %#ok<NASGU>
            loaded = TensorSpline.annotatedClassFromFile(path);
            tq = linspace(-1, 1, 31)';

            testCase.verifyTrue(ncfile.hasGroupWithName('knotAxes'))
            testCase.verifyEqual(loaded.knotPoints, spline.knotPoints)
            testCase.verifyEqual(loaded.xi, spline.xi)
            testCase.verifyEqual(loaded(tq), spline(tq), AbsTol=1e-12)
        end

        function tensorSplineTwoDimensionalRoundTripUsesKnotAxes(testCase)
            x = linspace(-1, 1, 6)';
            y = linspace(0, 2, 7)';
            knotPoints = {
                BSpline.knotPointsForDataPoints(x, S=3)
                BSpline.knotPointsForDataPoints(y, S=2)
            };
            xi = reshape(1:((numel(knotPoints{1}) - 4) * (numel(knotPoints{2}) - 3)), [], 1);
            spline = TensorSpline.fromKnotPoints(knotPoints, xi, S=[3 2], xMean=-2, xStd=3);

            path = string(strcat(tempname, '.nc'));
            cleanupPath = onCleanup(@() deleteIfPresent(path)); %#ok<NASGU>
            spline.writeToFile(path, shouldOverwriteExisting=true);

            loaded = TensorSpline.annotatedClassFromFile(path);
            [Xq, Yq] = ndgrid(linspace(-1, 1, 9)', linspace(0, 2, 11)');

            testCase.verifyEqual(loaded.knotPoints, spline.knotPoints)
            testCase.verifyEqual(loaded.xi, spline.xi)
            testCase.verifyEqual(loaded(Xq, Yq), spline(Xq, Yq), AbsTol=1e-12)
        end

        function interpolatingSplineOneDimensionalRoundTripPreservesGridVectors(testCase)
            t = linspace(0, 1, 11)';
            values = sin(2*pi*t);
            spline = InterpolatingSpline.fromGriddedValues(t, values, S=3);

            path = string(strcat(tempname, '.nc'));
            cleanupPath = onCleanup(@() deleteIfPresent(path)); %#ok<NASGU>
            spline.writeToFile(path, shouldOverwriteExisting=true);

            ncfile = NetCDFFile(path, shouldReadOnly=true);
            cleanupFile = onCleanup(@() ncfile.close()); %#ok<NASGU>
            loaded = InterpolatingSpline.annotatedClassFromFile(path);
            tq = linspace(0, 1, 51)';

            testCase.verifyTrue(ncfile.hasGroupWithName('gridAxes'))
            testCase.verifyFalse(any(strcmp({ncfile.realVariables.name}, 'values')))
            testCase.verifyEqual(loaded.gridVectors, spline.gridVectors)
            testCase.verifyEqual(loaded(tq), spline(tq), AbsTol=1e-12)
        end

        function interpolatingSplineTwoDimensionalRoundTripPreservesGridVectors(testCase)
            x = linspace(-1, 1, 6)';
            y = linspace(0, 2, 7)';
            [X, Y] = ndgrid(x, y);
            values = X.^2 .* Y.^3 + 2*X.*Y - 5;
            spline = InterpolatingSpline.fromGriddedValues({x, y}, values, S=[3 3]);

            path = string(strcat(tempname, '.nc'));
            cleanupPath = onCleanup(@() deleteIfPresent(path)); %#ok<NASGU>
            spline.writeToFile(path, shouldOverwriteExisting=true);

            loaded = InterpolatingSpline.annotatedClassFromFile(path);

            testCase.verifyEqual(loaded.gridVectors, spline.gridVectors)
            testCase.verifyEqual(loaded(X, Y), spline(X, Y), AbsTol=1e-10)
        end

        function constrainedSplineRoundTripRebuildsUnconstrainedDiagnostics(testCase)
            t = linspace(-1, 1, 21)';
            values = sin(pi*t) + 0.05*cos(3*pi*t);
            spline = ConstrainedSpline.fromData(t, values, S=3, distribution=NormalDistribution(sigma=1));
            originalSmoothing = spline.smoothingMatrix();

            path = string(strcat(tempname, '.nc'));
            cleanupPath = onCleanup(@() deleteIfPresent(path)); %#ok<NASGU>
            spline.writeToFile(path, shouldOverwriteExisting=true);

            ncfile = NetCDFFile(path, shouldReadOnly=true);
            cleanupFile = onCleanup(@() ncfile.close()); %#ok<NASGU>
            loaded = ConstrainedSpline.annotatedClassFromFile(path);
            tq = linspace(-1, 1, 41)';

            testCase.verifyTrue(ncfile.hasGroupWithName('gridAxes'))
            testCase.verifyTrue(ncfile.hasVariableWithName('dataPoints'))
            testCase.verifyTrue(ncfile.hasVariableWithName('dataValues'))
            testCase.verifyFalse(ncfile.hasVariableWithName('X'))
            testCase.verifyFalse(ncfile.hasVariableWithName('CmInv'))
            testCase.verifyFalse(ncfile.hasVariableWithName('W'))
            testCase.verifyFalse(ncfile.hasVariableWithName('Aeq'))
            testCase.verifyFalse(ncfile.hasVariableWithName('beq'))
            testCase.verifyFalse(ncfile.hasVariableWithName('Aineq'))
            testCase.verifyFalse(ncfile.hasVariableWithName('bineq'))
            testCase.verifyEqual(loaded.gridVectors, spline.gridVectors)
            testCase.verifyEqual(loaded.dataPoints, spline.dataPoints)
            testCase.verifyEqual(loaded.dataValues, spline.dataValues)
            testCase.verifyEqual(loaded(tq), spline(tq), AbsTol=1e-10)
            testCase.verifyEqual(loaded.X, spline.X, AbsTol=1e-12)
            testCase.verifyEqual(loaded.smoothingMatrix(), originalSmoothing, AbsTol=1e-10)
        end

        function constrainedSplineRoundTripPreservesConstraintsAndCorrelation(testCase)
            t = linspace(0, 1, 25)';
            values = 0.2 + exp(-20*(t - 0.5).^2);
            distribution = NormalDistribution(sigma=1);
            distribution.rho = @(tau) exp(-(tau/0.15).^2);
            constraints = [
                PointConstraint.equal(0.5, D=1, value=0)
                GlobalConstraint.positive()
            ];
            spline = ConstrainedSpline.fromData(t, values, S=3, distribution=distribution, constraints=constraints);

            path = string(strcat(tempname, '.nc'));
            cleanupPath = onCleanup(@() deleteIfPresent(path)); %#ok<NASGU>
            spline.writeToFile(path, shouldOverwriteExisting=true);

            loaded = ConstrainedSpline.annotatedClassFromFile(path);
            tq = linspace(0, 1, 51)';

            testCase.verifyEqual(func2str(loaded.distribution.rho), func2str(distribution.rho))
            testCase.verifyEqual(loaded.gridVectors, spline.gridVectors)
            testCase.verifyEqual(loaded.pointConstraints.points, spline.pointConstraints.points)
            testCase.verifyEqual(loaded.pointConstraints.D, spline.pointConstraints.D)
            testCase.verifyEqual(loaded.pointConstraints.value, spline.pointConstraints.value)
            testCase.verifyEqual(loaded.globalConstraints.shape, spline.globalConstraints.shape)
            testCase.verifyEqual(loaded(tq), spline(tq), AbsTol=1e-10)
            testCase.verifyError(@() loaded.smoothingMatrix(), 'ConstrainedSpline:UnavailableSmoothingMatrix')
        end

        function pointConstraintRoundTripPersistsScientificState(testCase)
            points = [0 0; 1 1; 2 2];
            constraint = PointConstraint.lowerBound(points, D=[1 0], value=[0; 1; 2]);

            path = string(strcat(tempname, '.nc'));
            cleanupPath = onCleanup(@() deleteIfPresent(path)); %#ok<NASGU>
            constraint.writeToFile(path, shouldOverwriteExisting=true);

            loaded = PointConstraint.annotatedClassFromFile(path);

            testCase.verifyEqual(loaded.points, constraint.points)
            testCase.verifyEqual(loaded.D, constraint.D)
            testCase.verifyEqual(loaded.relation, constraint.relation)
            testCase.verifyEqual(loaded.value, constraint.value)
        end

        function globalConstraintRoundTripPersistsScientificState(testCase)
            positiveConstraint = GlobalConstraint.positive();
            monotonicConstraint = GlobalConstraint.monotonicDecreasing(dimension=2);

            positivePath = string(strcat(tempname, '.nc'));
            cleanupPositivePath = onCleanup(@() deleteIfPresent(positivePath)); %#ok<NASGU>
            positiveConstraint.writeToFile(positivePath, shouldOverwriteExisting=true);
            loadedPositive = GlobalConstraint.annotatedClassFromFile(positivePath);

            monotonicPath = string(strcat(tempname, '.nc'));
            cleanupMonotonicPath = onCleanup(@() deleteIfPresent(monotonicPath)); %#ok<NASGU>
            monotonicConstraint.writeToFile(monotonicPath, shouldOverwriteExisting=true);
            loadedMonotonic = GlobalConstraint.annotatedClassFromFile(monotonicPath);

            testCase.verifyEqual(loadedPositive.shape, positiveConstraint.shape)
            testCase.verifyEmpty(loadedPositive.dimension)
            testCase.verifyEqual(loadedMonotonic.shape, monotonicConstraint.shape)
            testCase.verifyEqual(loadedMonotonic.dimension, monotonicConstraint.dimension)
        end

        function trajectorySplineConstructorRoundTripPreservesComponentSplines(testCase)
            t = linspace(0, 1, 17)';
            xSpline = ConstrainedSpline.fromData(t, cos(2*pi*t), S=3);
            ySpline = ConstrainedSpline.fromData(t, sin(2*pi*t), S=3);
            trajectory = TrajectorySpline(t=t, x=xSpline, y=ySpline);

            path = string(strcat(tempname, '.nc'));
            cleanupPath = onCleanup(@() deleteIfPresent(path)); %#ok<NASGU>
            trajectory.writeToFile(path, shouldOverwriteExisting=true);

            ncfile = NetCDFFile(path, shouldReadOnly=true);
            cleanupFile = onCleanup(@() ncfile.close()); %#ok<NASGU>
            loaded = TrajectorySpline.annotatedClassFromFile(path);

            testCase.verifyTrue(ncfile.hasGroupWithName('x'))
            testCase.verifyTrue(ncfile.hasGroupWithName('y'))
            testCase.verifyEqual(loaded.t, trajectory.t)
            testCase.verifyClass(loaded.x, 'ConstrainedSpline')
            testCase.verifyClass(loaded.y, 'ConstrainedSpline')
            testCase.verifyEqual(loaded.x(t), trajectory.x(t), AbsTol=1e-10)
            testCase.verifyEqual(loaded.y(t), trajectory.y(t), AbsTol=1e-10)
            testCase.verifyEqual(loaded.u(t), trajectory.u(t), AbsTol=1e-10)
            testCase.verifyEqual(loaded.v(t), trajectory.v(t), AbsTol=1e-10)
        end

        function trajectorySplineFromDataRoundTripPreservesSolvedState(testCase)
            t = linspace(0, 1, 17)';
            x = cos(2*pi*t);
            y = sin(2*pi*t);
            trajectory = TrajectorySpline.fromData(t, x, y, S=3);

            path = string(strcat(tempname, '.nc'));
            cleanupPath = onCleanup(@() deleteIfPresent(path)); %#ok<NASGU>
            trajectory.writeToFile(path, shouldOverwriteExisting=true);

            loaded = TrajectorySpline.annotatedClassFromFile(path);

            testCase.verifyEqual(loaded.t, trajectory.t)
            testCase.verifyClass(loaded.x, 'ConstrainedSpline')
            testCase.verifyClass(loaded.y, 'ConstrainedSpline')
            testCase.verifyEqual(loaded.x(t), trajectory.x(t), AbsTol=1e-10)
            testCase.verifyEqual(loaded.y(t), trajectory.y(t), AbsTol=1e-10)
            testCase.verifyEqual(loaded.u(t), trajectory.u(t), AbsTol=1e-10)
            testCase.verifyEqual(loaded.v(t), trajectory.v(t), AbsTol=1e-10)
        end
    end
end

function deleteIfPresent(path)
if isfile(path)
    delete(path);
end
end
