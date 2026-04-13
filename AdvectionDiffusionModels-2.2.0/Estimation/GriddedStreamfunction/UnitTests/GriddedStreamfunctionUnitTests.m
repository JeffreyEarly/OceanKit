classdef GriddedStreamfunctionUnitTests < matlab.unittest.TestCase

    methods (Test)
        function constructorFitsAndReturnsHandle(testCase)
            [~, t, x, y, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();

            fit = GriddedStreamfunction.fromTrajectories(trajectories);

            testCase.verifyClass(fit, "GriddedStreamfunction")
            testCase.verifyTrue(isa(fit, "handle"))
            testCase.verifyNotEmpty(fit.streamfunctionSpline)
            testCase.verifyNotEmpty(fit.observedTrajectories)
            testCase.verifyNotEmpty(fit.centerOfMassTrajectory)
            testCase.verifyNotEmpty(fit.backgroundTrajectory)
            testCase.verifyTrue(isstruct(fit.decomposition))
            testCase.verifyTrue(isfield(fit.decomposition, "fixedFrame"))
            testCase.verifyTrue(isfield(fit.decomposition, "centeredFrame"))
            testCase.verifyEqual(numel(fit.observedTrajectories), numel(trajectories))
            testCase.verifyEqual(numel(fit.decomposition.fixedFrame.background), numel(trajectories))
            testCase.verifyEqual(numel(fit.decomposition.fixedFrame.mesoscale), numel(trajectories))
            testCase.verifyEqual(numel(fit.decomposition.fixedFrame.submesoscale), numel(trajectories))
            testCase.verifyEqual(numel(fit.decomposition.centeredFrame.mesoscale), numel(trajectories))
            testCase.verifyEqual(numel(fit.decomposition.centeredFrame.submesoscale), numel(trajectories))
            testCase.verifyEqual(fit.fastS, 3)
            testCase.verifyEqual(fit.mesoscaleConstraint, "none")
            testCase.verifyGreaterThanOrEqual(fit.mesoscaleDegreesOfFreedom, 0)
            testCase.verifyLessThanOrEqual(max(abs(fit.centerOfMassTrajectory.x(t) - mean(x, 2))), 2e-1)
            testCase.verifyLessThanOrEqual(max(abs(fit.centerOfMassTrajectory.y(t) - mean(y, 2))), 2e-1)

            for iTrajectory = 1:numel(trajectories)
                testCase.verifyEqual(fit.observedTrajectories(iTrajectory).t, trajectories(iTrajectory).t)
            end
        end

        function removedPropertiesAreNotStored(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();

            fit = GriddedStreamfunction.fromTrajectories(trajectories);
            removedProperties = [ ...
                "sampleCounts", "observationTimes", "observedX", "observedY", ...
                "observedXVelocity", "observedYVelocity", "centeredX", "centeredY", ...
                "centerVelocityX", "centerVelocityY", "uMesoscaleObserved", ...
                "vMesoscaleObserved", "uMesoscaleComObserved", "vMesoscaleComObserved", ...
                "uMesoscaleRelativeObserved", "vMesoscaleRelativeObserved", ...
                "uBackgroundObserved", "vBackgroundObserved", "uSubmesoscaleObserved", ...
                "vSubmesoscaleObserved", "backgroundTrajectories", "mesoscaleTrajectories", ...
                "submesoscaleTrajectories", "centeredMesoscaleTrajectories", ...
                "centeredSubmesoscaleTrajectories", "fitDiagnostics"];

            for iProperty = 1:numel(removedProperties)
                testCase.verifyFalse(isprop(fit, removedProperties(iProperty)))
            end
        end

        function persistenceRoundTripPreservesSolvedState(testCase)
            [~, t, x, y, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            fit = GriddedStreamfunction.fromTrajectories(trajectories, mesoscaleConstraint="zeroVorticity");
            path = [tempname, '.nc'];
            cleanup = onCleanup(@() GriddedStreamfunctionUnitTests.deleteIfExists(path)); %#ok<NASGU>

            fit.writeToFile(path, shouldOverwriteExisting=true);
            restored = GriddedStreamfunction.fromFile(path);
            tGrid = repmat(t, 1, size(x, 2));

            testCase.verifyEqual(restored.mesoscaleConstraint, fit.mesoscaleConstraint)
            testCase.verifyEqual(restored.fastKnotPoints, fit.fastKnotPoints, "AbsTol", 1e-12)
            testCase.verifyEqual(restored.fastS, fit.fastS, "AbsTol", 1e-12)
            testCase.verifyEqual(restored.psiS, fit.psiS, "AbsTol", 1e-12)
            testCase.verifyEqual(restored.representativeTimes, fit.representativeTimes, "AbsTol", 1e-12)
            testCase.verifyEqual(restored.fitSupportTimes, fit.fitSupportTimes, "AbsTol", 1e-12)
            testCase.verifyEqual(restored.mesoscaleDegreesOfFreedom, fit.mesoscaleDegreesOfFreedom)
            for iDim = 1:3
                testCase.verifyEqual(restored.psiKnotPoints{iDim}, fit.psiKnotPoints{iDim}, "AbsTol", 1e-12)
            end

            testCase.verifyEqual(restored.backgroundTrajectory.x(t), fit.backgroundTrajectory.x(t), "AbsTol", 1e-12)
            testCase.verifyEqual(restored.backgroundTrajectory.y(t), fit.backgroundTrajectory.y(t), "AbsTol", 1e-12)
            testCase.verifyEqual(restored.uMesoscale(tGrid, x, y), fit.uMesoscale(tGrid, x, y), "AbsTol", 1e-12)
            testCase.verifyEqual(restored.vMesoscale(tGrid, x, y), fit.vMesoscale(tGrid, x, y), "AbsTol", 1e-12)
            testCase.verifyEqual(restored.sigma_n(tGrid, x, y), fit.sigma_n(tGrid, x, y), "AbsTol", 1e-12)
            testCase.verifyEqual(restored.sigma_s(tGrid, x, y), fit.sigma_s(tGrid, x, y), "AbsTol", 1e-12)
            testCase.verifyEqual(restored.zeta(tGrid, x, y), fit.zeta(tGrid, x, y), "AbsTol", 1e-12)

            for iTrajectory = 1:numel(trajectories)
                ti = fit.observedTrajectories(iTrajectory).t;
                GriddedStreamfunctionUnitTests.verifyTrajectorySplineEqualAtSamples( ...
                    testCase, restored.decomposition.fixedFrame.background(iTrajectory), ...
                    fit.decomposition.fixedFrame.background(iTrajectory), ti)
                GriddedStreamfunctionUnitTests.verifyTrajectorySplineEqualAtSamples( ...
                    testCase, restored.decomposition.fixedFrame.mesoscale(iTrajectory), ...
                    fit.decomposition.fixedFrame.mesoscale(iTrajectory), ti)
                GriddedStreamfunctionUnitTests.verifyTrajectorySplineEqualAtSamples( ...
                    testCase, restored.decomposition.fixedFrame.submesoscale(iTrajectory), ...
                    fit.decomposition.fixedFrame.submesoscale(iTrajectory), ti)
            end
        end

        function constructorRejectsRawSampleInputs(testCase)
            x = randn(10, 3);

            testCase.verifyError(@() GriddedStreamfunction.fromTrajectories(x), "MATLAB:validators:mustBeA")
        end

        function constructorRejectsCellArrayInputs(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            cellInput = cell(numel(trajectories), 1);
            for iTrajectory = 1:numel(trajectories)
                cellInput{iTrajectory} = trajectories(iTrajectory);
            end

            testCase.verifyError(@() GriddedStreamfunction.fromTrajectories(cellInput), "MATLAB:validators:mustBeA")
        end

        function constructorStoresMesoscaleConstraint(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();

            fit = GriddedStreamfunction.fromTrajectories(trajectories, mesoscaleConstraint="zeroVorticity");

            testCase.verifyEqual(fit.mesoscaleConstraint, "zeroVorticity")
        end

        function constructorRejectsInvalidMesoscaleConstraint(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();

            testCase.verifyError(@() GriddedStreamfunction.fromTrajectories(trajectories, mesoscaleConstraint="badConstraint"), "MATLAB:validators:mustBeMember")
        end

        function explicitNoneConstraintMatchesDefaultFit(testCase)
            [~, t, x, y, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();

            fitDefault = GriddedStreamfunction.fromTrajectories(trajectories);
            fitNone = GriddedStreamfunction.fromTrajectories(trajectories, mesoscaleConstraint="none");
            tGrid = repmat(t, 1, size(x, 2));

            testCase.verifyEqual(fitNone.mesoscaleConstraint, "none")
            testCase.verifyEqual(fitNone.streamfunctionSpline.xi, fitDefault.streamfunctionSpline.xi, "AbsTol", 1e-12)
            testCase.verifyEqual(fitNone.backgroundTrajectory.x(t), fitDefault.backgroundTrajectory.x(t), "AbsTol", 1e-12)
            testCase.verifyEqual(fitNone.backgroundTrajectory.y(t), fitDefault.backgroundTrajectory.y(t), "AbsTol", 1e-12)
            testCase.verifyEqual(fitNone.uMesoscale(tGrid, x, y), fitDefault.uMesoscale(tGrid, x, y), "AbsTol", 1e-12)
            testCase.verifyEqual(fitNone.vMesoscale(tGrid, x, y), fitDefault.vMesoscale(tGrid, x, y), "AbsTol", 1e-12)
        end

        function unconstrainedMesoscaleDegreesOfFreedomMatchesReducedBasisSize(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            fit = GriddedStreamfunction.fromTrajectories(trajectories, mesoscaleConstraint="none");
            [expectedDof, G, Aeq] = GriddedStreamfunctionUnitTests.expectedMesoscaleDegreesOfFreedom(fit);

            testCase.verifyEqual(size(Aeq, 1), 0)
            testCase.verifyEqual(fit.mesoscaleDegreesOfFreedom, expectedDof)
            testCase.verifyEqual(expectedDof, size(G, 2))
        end

        function zeroVorticityMesoscaleDegreesOfFreedomMatchesConstrainedCount(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            fit = GriddedStreamfunction.fromTrajectories(trajectories, mesoscaleConstraint="zeroVorticity");
            [expectedDof, G, Aeq] = GriddedStreamfunctionUnitTests.expectedMesoscaleDegreesOfFreedom(fit);

            testCase.verifyGreaterThan(size(Aeq, 1), 0)
            testCase.verifyEqual(fit.mesoscaleDegreesOfFreedom, expectedDof)
            testCase.verifyEqual(expectedDof, size(G, 2) - size(Aeq, 1))
        end

        function zeroStrainMesoscaleDegreesOfFreedomMatchesConstrainedCount(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            fit = GriddedStreamfunction.fromTrajectories(trajectories, mesoscaleConstraint="zeroStrain");
            [expectedDof, G, Aeq] = GriddedStreamfunctionUnitTests.expectedMesoscaleDegreesOfFreedom(fit);

            testCase.verifyGreaterThan(size(Aeq, 1), 0)
            testCase.verifyEqual(fit.mesoscaleDegreesOfFreedom, expectedDof)
            testCase.verifyEqual(expectedDof, size(G, 2) - size(Aeq, 1))
        end

        function visualPrincipalStrainAnglePreservesPrincipalValues(testCase)
            thetaExpected = [-30; -5; 20; 35];
            sigmaN = cosd(2 * thetaExpected);
            sigmaS = sind(2 * thetaExpected);

            thetaActual = GriddedStreamfunction.visualPrincipalStrainAngle(sigmaN, sigmaS);

            testCase.verifyEqual(thetaActual, thetaExpected, "AbsTol", 1e-12)
        end

        function visualPrincipalStrainAngleAvoidsBoundaryJump(testCase)
            thetaExpected = [35; 40; 45; 50; 55];
            sigmaN = cosd(2 * thetaExpected);
            sigmaS = sind(2 * thetaExpected);

            thetaActual = GriddedStreamfunction.visualPrincipalStrainAngle(sigmaN, sigmaS);

            testCase.verifyEqual(thetaActual, thetaExpected, "AbsTol", 1e-12)
        end

        function visualPrincipalStrainAngleReentersPrincipalRange(testCase)
            thetaExpected = [35; 50; 70; 40; 20];
            sigmaN = cosd(2 * thetaExpected);
            sigmaS = sind(2 * thetaExpected);

            thetaActual = GriddedStreamfunction.visualPrincipalStrainAngle(sigmaN, sigmaS);

            testCase.verifyEqual(thetaActual, thetaExpected, "AbsTol", 1e-12)
        end

        function visualPrincipalStrainAngleProcessesMatrixColumnsIndependently(testCase)
            thetaExpected = [40 -40; 50 -50; 60 -60; 30 -30];
            sigmaN = cosd(2 * thetaExpected);
            sigmaS = sind(2 * thetaExpected);

            thetaActual = GriddedStreamfunction.visualPrincipalStrainAngle(sigmaN, sigmaS);

            testCase.verifyEqual(thetaActual, thetaExpected, "AbsTol", 1e-12)
        end

        function visualPrincipalStrainAngleRestartsAfterNaNGap(testCase)
            thetaExpected = [35; 50; NaN; -40; -35];
            sigmaN = cosd(2 * thetaExpected);
            sigmaS = sind(2 * thetaExpected);
            sigmaN(3) = NaN;
            sigmaS(3) = NaN;

            thetaActual = GriddedStreamfunction.visualPrincipalStrainAngle(sigmaN, sigmaS);

            testCase.verifyEqual(thetaActual([1 2 4 5]), thetaExpected([1 2 4 5]), "AbsTol", 1e-12)
            testCase.verifyTrue(isnan(thetaActual(3)))
        end

        function visualPrincipalStrainAngleSupportsRadians(testCase)
            thetaExpected = [0.6; 0.9; 1.1];
            sigmaN = cos(2 * thetaExpected);
            sigmaS = sin(2 * thetaExpected);

            thetaActual = GriddedStreamfunction.visualPrincipalStrainAngle(sigmaN, sigmaS, units="radians");

            testCase.verifyEqual(thetaActual, thetaExpected, "AbsTol", 1e-12)
        end

        function synchronousLinearFieldRecovery(testCase)
            [model, t, x, y, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();

            fit = GriddedStreamfunction.fromTrajectories(trajectories);
            tGrid = repmat(t, 1, size(x, 2));
            uBackgroundGrid = repmat(fit.uBackground(t), 1, size(x, 2));
            vBackgroundGrid = repmat(fit.vBackground(t), 1, size(x, 2));
            uResolved = fit.uMesoscale(tGrid, x, y) + uBackgroundGrid;
            vResolved = fit.vMesoscale(tGrid, x, y) + vBackgroundGrid;

            testCase.verifyEqual(fit.representativeTimes, t, "AbsTol", 1e-12)
            testCase.verifyLessThanOrEqual(max(abs(fit.centerOfMassTrajectory.x(t) - mean(x, 2))), 2e-1)
            testCase.verifyLessThanOrEqual(max(abs(fit.centerOfMassTrajectory.y(t) - mean(y, 2))), 2e-1)
            testCase.verifyLessThanOrEqual(max(abs(fit.sigma_n(tGrid, x, y) - model.sigma_n), [], "all"), 4e-7)
            testCase.verifyLessThanOrEqual(max(abs(fit.sigma_s(tGrid, x, y) - model.sigma_s), [], "all"), 3e-7)
            testCase.verifyLessThanOrEqual(max(abs(fit.zeta(tGrid, x, y) - model.zeta), [], "all"), 2e-7)
            testCase.verifyLessThanOrEqual(max(abs(uResolved - model.u(tGrid, x, y)), [], "all"), 3e-4)
            testCase.verifyLessThanOrEqual(max(abs(vResolved - model.v(tGrid, x, y)), [], "all"), 3e-4)
        end

        function asynchronousLinearFieldRecovery(testCase)
            [model, trajectories] = GriddedStreamfunctionUnitTests.asynchronousLinearFieldData();
            fit = GriddedStreamfunction.fromTrajectories(trajectories);

            uResolvedError = zeros(numel(trajectories), 1);
            vResolvedError = zeros(numel(trajectories), 1);
            sigmaNError = zeros(numel(trajectories), 1);
            sigmaSError = zeros(numel(trajectories), 1);
            zetaError = zeros(numel(trajectories), 1);

            for iTrajectory = 1:numel(trajectories)
                trajectory = trajectories(iTrajectory);
                ti = trajectory.t;
                xi = trajectory.x(ti);
                yi = trajectory.y(ti);

                uResolvedError(iTrajectory) = max(abs(fit.uMesoscale(ti, xi, yi) + fit.uBackground(ti) - model.u(ti, xi, yi)));
                vResolvedError(iTrajectory) = max(abs(fit.vMesoscale(ti, xi, yi) + fit.vBackground(ti) - model.v(ti, xi, yi)));
                sigmaNError(iTrajectory) = max(abs(fit.sigma_n(ti, xi, yi) - model.sigma_n));
                sigmaSError(iTrajectory) = max(abs(fit.sigma_s(ti, xi, yi) - model.sigma_s));
                zetaError(iTrajectory) = max(abs(fit.zeta(ti, xi, yi) - model.zeta));
            end

            testCase.verifyLessThanOrEqual(max(uResolvedError), 5e-4)
            testCase.verifyLessThanOrEqual(max(vResolvedError), 5e-4)
            testCase.verifyLessThanOrEqual(max(sigmaNError), 2e-7)
            testCase.verifyLessThanOrEqual(max(sigmaSError), 2e-7)
            testCase.verifyLessThanOrEqual(max(zetaError), 2e-7)
        end

        function zeroVorticityConstraintZerosOffSampleVorticity(testCase)
            [~, trajectories, tQuery, xQuery, yQuery] = GriddedStreamfunctionUnitTests.offSampleLinearFieldData(4.0e-6, pi/9, 0);

            fit = GriddedStreamfunction.fromTrajectories(trajectories, mesoscaleConstraint="zeroVorticity");
            zetaValues = fit.zeta(tQuery, xQuery, yQuery);

            testCase.verifyEqual(fit.mesoscaleConstraint, "zeroVorticity")
            testCase.verifyTrue(all(isfinite(zetaValues), "all"))
            testCase.verifyLessThanOrEqual(max(abs(zetaValues), [], "all"), 1e-10)
            GriddedStreamfunctionUnitTests.verifyFiniteReconstruction(testCase, fit)
        end

        function zeroStrainConstraintZerosOffSampleStrain(testCase)
            [~, trajectories, tQuery, xQuery, yQuery] = GriddedStreamfunctionUnitTests.offSampleLinearFieldData(0, pi/9, -1.5e-6);

            fit = GriddedStreamfunction.fromTrajectories(trajectories, mesoscaleConstraint="zeroStrain");
            sigmaNValues = fit.sigma_n(tQuery, xQuery, yQuery);
            sigmaSValues = fit.sigma_s(tQuery, xQuery, yQuery);

            testCase.verifyEqual(fit.mesoscaleConstraint, "zeroStrain")
            testCase.verifyTrue(all(isfinite(sigmaNValues), "all"))
            testCase.verifyTrue(all(isfinite(sigmaSValues), "all"))
            testCase.verifyLessThanOrEqual(max(abs(sigmaNValues), [], "all"), 1e-10)
            testCase.verifyLessThanOrEqual(max(abs(sigmaSValues), [], "all"), 1e-10)
            GriddedStreamfunctionUnitTests.verifyFiniteReconstruction(testCase, fit)
        end

        function anchoredTrajectoryComponentsMatchConventionsAndReconstructPositions(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            fit = GriddedStreamfunction.fromTrajectories(trajectories);

            for iTrajectory = 1:numel(trajectories)
                trajectory = fit.observedTrajectories(iTrajectory);
                ti = trajectory.t;
                observedX = trajectory.x(ti);
                observedY = trajectory.y(ti);
                [~, centeredX, centeredY] = fit.centeredCoordinates(ti, observedX, observedY);

                backgroundX = fit.decomposition.fixedFrame.background(iTrajectory).x(ti);
                backgroundY = fit.decomposition.fixedFrame.background(iTrajectory).y(ti);
                mesoscaleX = fit.decomposition.fixedFrame.mesoscale(iTrajectory).x(ti);
                mesoscaleY = fit.decomposition.fixedFrame.mesoscale(iTrajectory).y(ti);
                submesoscaleX = fit.decomposition.fixedFrame.submesoscale(iTrajectory).x(ti);
                submesoscaleY = fit.decomposition.fixedFrame.submesoscale(iTrajectory).y(ti);
                centeredMesoscaleX = fit.decomposition.centeredFrame.mesoscale(iTrajectory).x(ti);
                centeredMesoscaleY = fit.decomposition.centeredFrame.mesoscale(iTrajectory).y(ti);
                centeredSubmesoscaleX = fit.decomposition.centeredFrame.submesoscale(iTrajectory).x(ti);
                centeredSubmesoscaleY = fit.decomposition.centeredFrame.submesoscale(iTrajectory).y(ti);
                backgroundShiftedX = fit.backgroundTrajectory.x(ti) - fit.backgroundTrajectory.x(ti(1));
                backgroundShiftedY = fit.backgroundTrajectory.y(ti) - fit.backgroundTrajectory.y(ti(1));

                testCase.verifyEqual(backgroundX(1), 0, "AbsTol", 1e-12)
                testCase.verifyEqual(backgroundY(1), 0, "AbsTol", 1e-12)
                testCase.verifyEqual(submesoscaleX(1), 0, "AbsTol", 1e-12)
                testCase.verifyEqual(submesoscaleY(1), 0, "AbsTol", 1e-12)
                testCase.verifyEqual(centeredSubmesoscaleX(1), 0, "AbsTol", 1e-12)
                testCase.verifyEqual(centeredSubmesoscaleY(1), 0, "AbsTol", 1e-12)
                testCase.verifyEqual(mesoscaleX(1), observedX(1), "AbsTol", 1e-12)
                testCase.verifyEqual(mesoscaleY(1), observedY(1), "AbsTol", 1e-12)
                testCase.verifyEqual(centeredMesoscaleX(1), centeredX(1), "AbsTol", 1e-12)
                testCase.verifyEqual(centeredMesoscaleY(1), centeredY(1), "AbsTol", 1e-12)
                testCase.verifyEqual(backgroundX, backgroundShiftedX, "AbsTol", 1e-12)
                testCase.verifyEqual(backgroundY, backgroundShiftedY, "AbsTol", 1e-12)

                testCase.verifyEqual(backgroundX + mesoscaleX + submesoscaleX, observedX, "AbsTol", 1e-2)
                testCase.verifyEqual(backgroundY + mesoscaleY + submesoscaleY, observedY, "AbsTol", 1e-2)
                testCase.verifyEqual(centeredMesoscaleX + centeredSubmesoscaleX, centeredX, "AbsTol", 1e-2)
                testCase.verifyEqual(centeredMesoscaleY + centeredSubmesoscaleY, centeredY, "AbsTol", 1e-2)
            end
        end

        function componentTrajectoryDerivativesMatchVelocityDecomposition(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            fit = GriddedStreamfunction.fromTrajectories(trajectories);

            for iTrajectory = 1:numel(trajectories)
                trajectory = fit.observedTrajectories(iTrajectory);
                ti = trajectory.t;
                observedX = trajectory.x(ti);
                observedY = trajectory.y(ti);
                observedU = trajectory.u(ti);
                observedV = trajectory.v(ti);

                background = fit.decomposition.fixedFrame.background(iTrajectory);
                mesoscale = fit.decomposition.fixedFrame.mesoscale(iTrajectory);
                submesoscale = fit.decomposition.fixedFrame.submesoscale(iTrajectory);
                centeredMesoscale = fit.decomposition.centeredFrame.mesoscale(iTrajectory);
                centeredSubmesoscale = fit.decomposition.centeredFrame.submesoscale(iTrajectory);

                testCase.verifyEqual(background.u(ti), fit.uBackground(ti), "AbsTol", 1e-12)
                testCase.verifyEqual(background.v(ti), fit.vBackground(ti), "AbsTol", 1e-12)
                testCase.verifyEqual(mesoscale.u(ti), fit.uMesoscale(ti, observedX, observedY), "AbsTol", 1e-12)
                testCase.verifyEqual(mesoscale.v(ti), fit.vMesoscale(ti, observedX, observedY), "AbsTol", 1e-12)
                testCase.verifyEqual(background.u(ti) + mesoscale.u(ti) + submesoscale.u(ti), observedU, "AbsTol", 1e-12)
                testCase.verifyEqual(background.v(ti) + mesoscale.v(ti) + submesoscale.v(ti), observedV, "AbsTol", 1e-12)
                testCase.verifyEqual( ...
                    centeredMesoscale.u(ti) + centeredSubmesoscale.u(ti), ...
                    observedU - fit.centerOfMassTrajectory.u(ti), ...
                    "AbsTol", 1e-12)
                testCase.verifyEqual( ...
                    centeredMesoscale.v(ti) + centeredSubmesoscale.v(ti), ...
                    observedV - fit.centerOfMassTrajectory.v(ti), ...
                    "AbsTol", 1e-12)
            end
        end

        function decomposeTrajectoriesReproducesStoredDecomposition(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            fit = GriddedStreamfunction.fromTrajectories(trajectories);
            decomposition = fit.decomposeTrajectories(fit.observedTrajectories);

            for iTrajectory = 1:numel(trajectories)
                ti = fit.observedTrajectories(iTrajectory).t;

                GriddedStreamfunctionUnitTests.verifyTrajectorySplineEqualAtSamples( ...
                    testCase, decomposition.fixedFrame.background(iTrajectory), ...
                    fit.decomposition.fixedFrame.background(iTrajectory), ti)
                GriddedStreamfunctionUnitTests.verifyTrajectorySplineEqualAtSamples( ...
                    testCase, decomposition.fixedFrame.mesoscale(iTrajectory), ...
                    fit.decomposition.fixedFrame.mesoscale(iTrajectory), ti)
                GriddedStreamfunctionUnitTests.verifyTrajectorySplineEqualAtSamples( ...
                    testCase, decomposition.fixedFrame.submesoscale(iTrajectory), ...
                    fit.decomposition.fixedFrame.submesoscale(iTrajectory), ti)
                GriddedStreamfunctionUnitTests.verifyTrajectorySplineEqualAtSamples( ...
                    testCase, decomposition.centeredFrame.mesoscale(iTrajectory), ...
                    fit.decomposition.centeredFrame.mesoscale(iTrajectory), ti)
                GriddedStreamfunctionUnitTests.verifyTrajectorySplineEqualAtSamples( ...
                    testCase, decomposition.centeredFrame.submesoscale(iTrajectory), ...
                    fit.decomposition.centeredFrame.submesoscale(iTrajectory), ti)
            end
        end

        function decomposeTrajectoriesSupportsTrimmedAsynchronousInputs(testCase)
            [fitTrajectories, evaluationTrajectories] = GriddedStreamfunctionUnitTests.trimmedTrajectoryData();
            fit = GriddedStreamfunction.fromTrajectories(fitTrajectories);
            decomposition = fit.decomposeTrajectories(evaluationTrajectories);

            testCase.verifyEqual(numel(decomposition.fixedFrame.background), numel(evaluationTrajectories))
            testCase.verifyEqual(numel(decomposition.fixedFrame.mesoscale), numel(evaluationTrajectories))
            testCase.verifyEqual(numel(decomposition.fixedFrame.submesoscale), numel(evaluationTrajectories))
            testCase.verifyEqual(numel(decomposition.centeredFrame.mesoscale), numel(evaluationTrajectories))
            testCase.verifyEqual(numel(decomposition.centeredFrame.submesoscale), numel(evaluationTrajectories))

            for iTrajectory = 1:numel(evaluationTrajectories)
                trajectory = evaluationTrajectories(iTrajectory);
                ti = trajectory.t;
                observedU = trajectory.u(ti);
                observedV = trajectory.v(ti);

                background = decomposition.fixedFrame.background(iTrajectory);
                mesoscale = decomposition.fixedFrame.mesoscale(iTrajectory);
                submesoscale = decomposition.fixedFrame.submesoscale(iTrajectory);
                centeredMesoscale = decomposition.centeredFrame.mesoscale(iTrajectory);
                centeredSubmesoscale = decomposition.centeredFrame.submesoscale(iTrajectory);

                testCase.verifyTrue(all(isfinite(background.x(ti))))
                testCase.verifyTrue(all(isfinite(background.y(ti))))
                testCase.verifyTrue(all(isfinite(mesoscale.x(ti))))
                testCase.verifyTrue(all(isfinite(mesoscale.y(ti))))
                testCase.verifyTrue(all(isfinite(submesoscale.x(ti))))
                testCase.verifyTrue(all(isfinite(submesoscale.y(ti))))
                testCase.verifyEqual(background.u(ti) + mesoscale.u(ti) + submesoscale.u(ti), observedU, "AbsTol", 1e-12)
                testCase.verifyEqual(background.v(ti) + mesoscale.v(ti) + submesoscale.v(ti), observedV, "AbsTol", 1e-12)
                testCase.verifyEqual( ...
                    centeredMesoscale.u(ti) + centeredSubmesoscale.u(ti), ...
                    observedU - fit.centerOfMassTrajectory.u(ti), ...
                    "AbsTol", 1e-12)
                testCase.verifyEqual( ...
                    centeredMesoscale.v(ti) + centeredSubmesoscale.v(ti), ...
                    observedV - fit.centerOfMassTrajectory.v(ti), ...
                    "AbsTol", 1e-12)
            end
        end

        function decomposeTrajectoriesRejectsTrajectoriesOutsideFitDomains(testCase)
            [~, t, x, y, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            fit = GriddedStreamfunction.fromTrajectories(trajectories);

            timeOutsideTrajectories = GriddedStreamfunctionUnitTests.trajectorySplinesFromMatrices(t + 2 * 900, x, y);
            qRange = fit.psiKnotPoints{1}(end) - fit.psiKnotPoints{1}(1);
            psiOutsideTrajectories = GriddedStreamfunctionUnitTests.trajectorySplinesFromMatrices(t, x + 2 * qRange, y);

            testCase.verifyError(@() fit.decomposeTrajectories(timeOutsideTrajectories), ...
                "GriddedStreamfunction:ObservationOutsideFastDomain")
            testCase.verifyError(@() fit.decomposeTrajectories(psiOutsideTrajectories), ...
                "GriddedStreamfunction:ObservationOutsidePsiDomain")
        end

        function constructorReadsTrajectoryEvaluationsInsteadOfDataValues(testCase)
            trajectories = GriddedStreamfunctionUnitTests.smoothedTrajectoryData();
            fit = GriddedStreamfunction.fromTrajectories(trajectories, psiKnotPoints=GriddedStreamfunctionUnitTests.zeroPsiKnotPoints());

            maxObservedRawDifference = 0;
            for iTrajectory = 1:numel(trajectories)
                trajectory = fit.observedTrajectories(iTrajectory);
                ti = trajectory.t;
                observedX = trajectory.x(ti);
                observedY = trajectory.y(ti);
                rawX = trajectory.x.dataValues;

                reconstructedX = fit.decomposition.fixedFrame.background(iTrajectory).x(ti) + ...
                    fit.decomposition.fixedFrame.mesoscale(iTrajectory).x(ti) + ...
                    fit.decomposition.fixedFrame.submesoscale(iTrajectory).x(ti);
                reconstructedY = fit.decomposition.fixedFrame.background(iTrajectory).y(ti) + ...
                    fit.decomposition.fixedFrame.mesoscale(iTrajectory).y(ti) + ...
                    fit.decomposition.fixedFrame.submesoscale(iTrajectory).y(ti);

                testCase.verifyEqual(reconstructedX, observedX, "AbsTol", 1e-2)
                testCase.verifyEqual(reconstructedY, observedY, "AbsTol", 1e-2)
                maxObservedRawDifference = max(maxObservedRawDifference, max(abs(observedX - rawX)));
            end

            testCase.verifyGreaterThan(maxObservedRawDifference, 1e-3)
        end

        function decompositionSampleDataMatchesObservedTrajectorySamples(testCase)
            trajectories = GriddedStreamfunctionUnitTests.smoothedTrajectoryData();
            fit = GriddedStreamfunction.fromTrajectories(trajectories, psiKnotPoints=GriddedStreamfunctionUnitTests.zeroPsiKnotPoints());
            sampleData = fit.decompositionSampleData(fit.observedTrajectories);
            expectedSamples = GriddedStreamfunctionUnitTests.directTrajectorySamples(fit.observedTrajectories);

            for iTrajectory = 1:numel(expectedSamples.tCell)
                testCase.verifyEqual(sampleData.tCell{iTrajectory}, expectedSamples.tCell{iTrajectory}, "AbsTol", 1e-12)
                testCase.verifyEqual(sampleData.xCell{iTrajectory}, expectedSamples.xCell{iTrajectory}, "AbsTol", 1e-12)
                testCase.verifyEqual(sampleData.yCell{iTrajectory}, expectedSamples.yCell{iTrajectory}, "AbsTol", 1e-12)
                testCase.verifyEqual(sampleData.xDotCell{iTrajectory}, expectedSamples.xDotCell{iTrajectory}, "AbsTol", 1e-12)
                testCase.verifyEqual(sampleData.yDotCell{iTrajectory}, expectedSamples.yDotCell{iTrajectory}, "AbsTol", 1e-12)
            end

            allT = expectedSamples.allT;
            allX = expectedSamples.allX;
            allY = expectedSamples.allY;
            allXDot = expectedSamples.allXDot;
            allYDot = expectedSamples.allYDot;
            B = reshape(BSpline.matrixForDataPoints(allT, knotPoints=fit.fastKnotPoints, S=fit.fastS), numel(allT), []);
            mxDotAll = reshape(fit.centerOfMassTrajectory.u(allT), [], 1);
            myDotAll = reshape(fit.centerOfMassTrajectory.v(allT), [], 1);
            mesoscaleObserved = [ ...
                reshape(fit.uMesoscale(allT, allX, allY), [], 1), ...
                reshape(fit.vMesoscale(allT, allX, allY), [], 1)];
            mesoscaleComObserved = B * (B \ mesoscaleObserved);
            uBackgroundObserved = mxDotAll - mesoscaleComObserved(:, 1);
            vBackgroundObserved = myDotAll - mesoscaleComObserved(:, 2);

            testCase.verifyEqual(sampleData.allT, allT, "AbsTol", 1e-12)
            testCase.verifyEqual(sampleData.allX, allX, "AbsTol", 1e-12)
            testCase.verifyEqual(sampleData.allY, allY, "AbsTol", 1e-12)
            testCase.verifyEqual(sampleData.allXDot, allXDot, "AbsTol", 1e-12)
            testCase.verifyEqual(sampleData.allYDot, allYDot, "AbsTol", 1e-12)
            testCase.verifyEqual(sampleData.sampleSupportTimes, expectedSamples.sampleSupportTimes, "AbsTol", 1e-12)
            testCase.verifyEqual(full(sampleData.B), full(B), "AbsTol", 1e-12)
            testCase.verifyEqual(sampleData.uMesoscaleComObserved, mesoscaleComObserved(:, 1), "AbsTol", 1e-12)
            testCase.verifyEqual(sampleData.vMesoscaleComObserved, mesoscaleComObserved(:, 2), "AbsTol", 1e-12)
            testCase.verifyEqual(sampleData.uSubmesoscaleObserved, allXDot - uBackgroundObserved - mesoscaleObserved(:, 1), "AbsTol", 1e-12)
            testCase.verifyEqual(sampleData.vSubmesoscaleObserved, allYDot - vBackgroundObserved - mesoscaleObserved(:, 2), "AbsTol", 1e-12)
        end

        function decompositionSampleDataUsesPassedTrajectoryOrdering(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            fit = GriddedStreamfunction.fromTrajectories(trajectories);
            reorderedTrajectories = flipud(fit.observedTrajectories);
            sampleData = fit.decompositionSampleData(reorderedTrajectories);
            expectedSamples = GriddedStreamfunctionUnitTests.directTrajectorySamples(reorderedTrajectories);

            for iTrajectory = 1:numel(expectedSamples.tCell)
                testCase.verifyEqual(sampleData.tCell{iTrajectory}, expectedSamples.tCell{iTrajectory}, "AbsTol", 1e-12)
                testCase.verifyEqual(sampleData.xCell{iTrajectory}, expectedSamples.xCell{iTrajectory}, "AbsTol", 1e-12)
                testCase.verifyEqual(sampleData.yCell{iTrajectory}, expectedSamples.yCell{iTrajectory}, "AbsTol", 1e-12)
                testCase.verifyEqual(sampleData.xDotCell{iTrajectory}, expectedSamples.xDotCell{iTrajectory}, "AbsTol", 1e-12)
                testCase.verifyEqual(sampleData.yDotCell{iTrajectory}, expectedSamples.yDotCell{iTrajectory}, "AbsTol", 1e-12)
            end

            testCase.verifyEqual(sampleData.allT, expectedSamples.allT, "AbsTol", 1e-12)
            testCase.verifyEqual(sampleData.allX, expectedSamples.allX, "AbsTol", 1e-12)
            testCase.verifyEqual(sampleData.allY, expectedSamples.allY, "AbsTol", 1e-12)
            testCase.verifyEqual(sampleData.allXDot, expectedSamples.allXDot, "AbsTol", 1e-12)
            testCase.verifyEqual(sampleData.allYDot, expectedSamples.allYDot, "AbsTol", 1e-12)
            testCase.verifyEqual(sampleData.sampleSupportTimes, expectedSamples.sampleSupportTimes, "AbsTol", 1e-12)
        end

        function decompositionSampleDataMesoscaleValuesMatchLegacyEvaluators(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            fit = GriddedStreamfunction.fromTrajectories(trajectories);
            sampleData = fit.decompositionSampleData(fit.observedTrajectories);

            expectedU = reshape(fit.uMesoscale(sampleData.allT, sampleData.allX, sampleData.allY), [], 1);
            expectedV = reshape(fit.vMesoscale(sampleData.allT, sampleData.allX, sampleData.allY), [], 1);

            testCase.verifyEqual(sampleData.uMesoscaleObserved, expectedU, "AbsTol", 1e-12)
            testCase.verifyEqual(sampleData.vMesoscaleObserved, expectedV, "AbsTol", 1e-12)
        end

        function resampledTrajectorySampleDataMatchesDirectSampling(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            observedSampleData = GriddedStreamfunction.sampleTrajectoryData(trajectories);
            sampledIndices = [3 1 3 2].';
            expectedSamples = GriddedStreamfunctionUnitTests.directTrajectorySamples(trajectories(sampledIndices));
            sampleData = GriddedStreamfunction.resampledTrajectorySampleData(observedSampleData, sampledIndices);

            for iTrajectory = 1:numel(expectedSamples.tCell)
                testCase.verifyEqual(sampleData.tCell{iTrajectory}, expectedSamples.tCell{iTrajectory}, "AbsTol", 1e-12)
                testCase.verifyEqual(sampleData.xCell{iTrajectory}, expectedSamples.xCell{iTrajectory}, "AbsTol", 1e-12)
                testCase.verifyEqual(sampleData.yCell{iTrajectory}, expectedSamples.yCell{iTrajectory}, "AbsTol", 1e-12)
                testCase.verifyEqual(sampleData.xDotCell{iTrajectory}, expectedSamples.xDotCell{iTrajectory}, "AbsTol", 1e-12)
                testCase.verifyEqual(sampleData.yDotCell{iTrajectory}, expectedSamples.yDotCell{iTrajectory}, "AbsTol", 1e-12)
            end

            testCase.verifyEqual(sampleData.allT, expectedSamples.allT, "AbsTol", 1e-12)
            testCase.verifyEqual(sampleData.allX, expectedSamples.allX, "AbsTol", 1e-12)
            testCase.verifyEqual(sampleData.allY, expectedSamples.allY, "AbsTol", 1e-12)
            testCase.verifyEqual(sampleData.allXDot, expectedSamples.allXDot, "AbsTol", 1e-12)
            testCase.verifyEqual(sampleData.allYDot, expectedSamples.allYDot, "AbsTol", 1e-12)
            testCase.verifyEqual(sampleData.sampleSupportTimes, expectedSamples.sampleSupportTimes, "AbsTol", 1e-12)
        end

        function reducedBasisGaugeRemovesScalarSpatialMode(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            fit = GriddedStreamfunction.fromTrajectories(trajectories);

            basisSize = fit.streamfunctionSpline.basisSize;
            xiByTime = reshape(fit.streamfunctionSpline.xi, prod(basisSize(1:2)), basisSize(3));
            gaugeModeMatrix = GriddedStreamfunctionUnitTests.constantGaugeMode(fit.psiKnotPoints, fit.psiS, basisSize);

            testCase.verifyEqual(size(gaugeModeMatrix, 2), 1)
            testCase.verifyEqual(gaugeModeMatrix.' * xiByTime, zeros(1, basisSize(3)), "AbsTol", 1e-10)
        end

        function fastBasisProjectionMatchesLegacyProjector(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            fit = GriddedStreamfunction.fromTrajectories(trajectories);

            nTrajectories = numel(trajectories);
            tCell = cell(nTrajectories, 1);
            xCell = cell(nTrajectories, 1);
            yCell = cell(nTrajectories, 1);
            for iTrajectory = 1:nTrajectories
                ti = reshape(trajectories(iTrajectory).t, [], 1);
                tCell{iTrajectory} = ti;
                xCell{iTrajectory} = reshape(trajectories(iTrajectory).x(ti), [], 1);
                yCell{iTrajectory} = reshape(trajectories(iTrajectory).y(ti), [], 1);
            end

            allT = vertcat(tCell{:});
            allX = vertcat(xCell{:});
            allY = vertcat(yCell{:});
            B = reshape(BSpline.matrixForDataPoints(allT, knotPoints=fit.fastKnotPoints, S=fit.fastS), numel(allT), []);
            dBTensor = BSpline.matrixForDataPoints(allT, knotPoints=fit.fastKnotPoints, S=fit.fastS, D=1);
            dB = reshape(dBTensor(:, :, 2), numel(allT), []);
            centerCoefficients = B \ [allX, allY];
            mxAll = B * centerCoefficients(:, 1);
            myAll = B * centerCoefficients(:, 2);
            mxDotAll = dB * centerCoefficients(:, 1);
            myDotAll = dB * centerCoefficients(:, 2);

            testCase.verifyEqual(reshape(fit.centerOfMassTrajectory.x(allT), [], 1), mxAll, "AbsTol", 1e-12)
            testCase.verifyEqual(reshape(fit.centerOfMassTrajectory.y(allT), [], 1), myAll, "AbsTol", 1e-12)
            testCase.verifyEqual(reshape(fit.centerOfMassTrajectory.u(allT), [], 1), mxDotAll, "AbsTol", 1e-12)
            testCase.verifyEqual(reshape(fit.centerOfMassTrajectory.v(allT), [], 1), myDotAll, "AbsTol", 1e-12)

            sampleData = fit.decompositionSampleData(trajectories);
            mesoscaleObserved = [ ...
                reshape(fit.uMesoscale(allT, allX, allY), [], 1), ...
                reshape(fit.vMesoscale(allT, allX, allY), [], 1)];
            mesoscaleComObserved = B * (B \ mesoscaleObserved);

            testCase.verifyEqual(full(sampleData.B), full(B), "AbsTol", 1e-12)
            testCase.verifyEqual(sampleData.uMesoscaleComObserved, mesoscaleComObserved(:, 1), "AbsTol", 1e-12)
            testCase.verifyEqual(sampleData.vMesoscaleComObserved, mesoscaleComObserved(:, 2), "AbsTol", 1e-12)
        end

        function evaluationInputShapes(testCase)
            [~, t, x, y, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            fit = GriddedStreamfunction.fromTrajectories(trajectories);

            tGrid = repmat(t, 1, size(x, 2));
            mxGrid = repmat(fit.centerOfMassTrajectory.x(t), 1, size(x, 2));
            myGrid = repmat(fit.centerOfMassTrajectory.y(t), 1, size(x, 2));
            scalarTime = t(4);
            scalarRow = repmat(scalarTime, 1, size(x, 2));

            [tEval, q, r] = fit.centeredCoordinates(t, x, y);

            testCase.verifyEqual(tEval, tGrid, "AbsTol", 1e-12)
            testCase.verifyEqual(q, x - mxGrid, "AbsTol", 1e-12)
            testCase.verifyEqual(r, y - myGrid, "AbsTol", 1e-12)
            testCase.verifyEqual(fit.uMesoscale(scalarTime, x(4, :), y(4, :)), fit.uMesoscale(scalarRow, x(4, :), y(4, :)), "AbsTol", 1e-12)
            testCase.verifyEqual(fit.psiMesoscale(t, x, y), fit.psiMesoscale(tGrid, x, y), "AbsTol", 1e-12)
            testCase.verifyEqual(fit.uMesoscale(t, x, y), fit.uMesoscale(tGrid, x, y), "AbsTol", 1e-12)
            testCase.verifyEqual(fit.vMesoscale(t, x, y), fit.vMesoscale(tGrid, x, y), "AbsTol", 1e-12)

            testCase.verifyError(@() fit.centeredCoordinates(t, x(1:(end - 1), :), y), "GriddedStreamfunction:EvaluationSizeMismatch")
            testCase.verifyError(@() fit.centeredCoordinates(t(1:(end - 1)), x, y), "GriddedStreamfunction:EvaluationTimeSizeMismatch")
            testCase.verifyError(@() fit.uMesoscale(t, x(1:(end - 1), :), y), "GriddedStreamfunction:EvaluationSizeMismatch")
            testCase.verifyError(@() fit.uMesoscale(t(1:(end - 1)), x, y), "GriddedStreamfunction:EvaluationTimeSizeMismatch")
        end
    end

    methods (Static, Access = private)
        function [model, t, x, y, trajectories] = synchronousLinearFieldData()
            sigma = 4.0e-6;
            theta = pi/9;
            zeta = -1.5e-6;
            u0 = 0.08;
            v0 = -0.05;

            model = LinearVelocityField(sigma=sigma, theta=theta, zeta=zeta, u0=u0, v0=v0);
            integrator = AdvectionDiffusionIntegrator(model, 0);

            [x0, y0] = ndgrid([-600 0 600] + 1200, [-400 0 400] - 800);
            [t, x, y] = integrator.particleTrajectories(x0(:), y0(:), 12*3600, 900);
            trajectories = GriddedStreamfunctionUnitTests.trajectorySplinesFromMatrices(t, x, y);
        end

        function [model, trajectories] = asynchronousLinearFieldData()
            sigma = 4.0e-6;
            theta = pi/9;
            zeta = -1.5e-6;
            u0 = 0.08;
            v0 = -0.05;

            model = LinearVelocityField(sigma=sigma, theta=theta, zeta=zeta, u0=u0, v0=v0);
            integrator = AdvectionDiffusionIntegrator(model, 0);

            [x0, y0] = ndgrid([-600 0 600] + 1200, [-400 0 400] - 800);
            [tFine, xFine, yFine] = integrator.particleTrajectories(x0(:), y0(:), 12*3600, 300);

            nDrifters = size(xFine, 2);
            trajectories = TrajectorySpline.empty(0, 1);
            drifterOffset = [-1; 0; 1; 1; 0; 1; 1; 0; -1];
            baseIndices = 3:3:(numel(tFine) - 2);
            jitterPattern = [0; 1; 0; -1; 0; 1; 0; -1; 0; 0];
            jitter = repmat(jitterPattern, ceil(numel(baseIndices) / numel(jitterPattern)), 1);
            jitter = jitter(1:numel(baseIndices));

            for iDrifter = 1:nDrifters
                indices = baseIndices + transpose(jitter) + drifterOffset(iDrifter);
                trajectories(end + 1, 1) = TrajectorySpline.fromData(tFine(indices), xFine(indices, iDrifter), yFine(indices, iDrifter), S=3); %#ok<AGROW>
            end
        end

        function [model, trajectories, tQuery, xQuery, yQuery] = offSampleLinearFieldData(sigma, theta, zeta)
            u0 = 0.08;
            v0 = -0.05;

            model = LinearVelocityField(sigma=sigma, theta=theta, zeta=zeta, u0=u0, v0=v0);
            integrator = AdvectionDiffusionIntegrator(model, 0);
            [x0, y0] = ndgrid([-600 0 600] + 1200, [-400 0 400] - 800);
            [tFine, xFine, yFine] = integrator.particleTrajectories(x0(:), y0(:), 12*3600, 300);

            coarseIndices = 1:3:numel(tFine);
            queryIndices = setdiff(1:numel(tFine), coarseIndices);
            trajectories = GriddedStreamfunctionUnitTests.trajectorySplinesFromMatrices( ...
                tFine(coarseIndices), xFine(coarseIndices, :), yFine(coarseIndices, :));

            tQuery = repmat(tFine(queryIndices), 1, size(xFine, 2));
            xQuery = xFine(queryIndices, :);
            yQuery = yFine(queryIndices, :);
        end

        function trajectories = trajectorySplinesFromMatrices(t, x, y, splineDegree)
            if nargin < 4
                splineDegree = 3;
            end

            trajectories = TrajectorySpline.empty(0, 1);
            for iDrifter = 1:size(x, 2)
                trajectories(end + 1, 1) = TrajectorySpline.fromData(t, x(:, iDrifter), y(:, iDrifter), S=splineDegree); %#ok<AGROW>
            end
        end

        function trajectories = smoothedTrajectoryData()
            t = linspace(0, 1, 21).';
            x1 = 2 + 0.3 * t;
            y1 = -1 + 0.2 * t;
            x2 = 3 + 0.3 * t;
            y2 = -0.5 + 0.2 * t;
            perturbation = 0.05 * sin(2*pi*t);

            xSpline1 = ConstrainedSpline.fromData(t, x1 + perturbation, S=3, splineDOF=6);
            ySpline1 = ConstrainedSpline.fromData(t, y1 - perturbation, S=3, splineDOF=6);
            xSpline2 = ConstrainedSpline.fromData(t, x2 - perturbation, S=3, splineDOF=6);
            ySpline2 = ConstrainedSpline.fromData(t, y2 + perturbation, S=3, splineDOF=6);

            trajectories = TrajectorySpline.empty(0, 1);
            trajectories(end + 1, 1) = TrajectorySpline(t=t, x=xSpline1, y=ySpline1);
            trajectories(end + 1, 1) = TrajectorySpline(t=t, x=xSpline2, y=ySpline2);
        end

        function sampleData = directTrajectorySamples(trajectories)
            nTrajectories = numel(trajectories);
            tCell = cell(nTrajectories, 1);
            xCell = cell(nTrajectories, 1);
            yCell = cell(nTrajectories, 1);
            xDotCell = cell(nTrajectories, 1);
            yDotCell = cell(nTrajectories, 1);

            for iTrajectory = 1:nTrajectories
                ti = reshape(trajectories(iTrajectory).t, [], 1);
                tCell{iTrajectory} = ti;
                xCell{iTrajectory} = reshape(trajectories(iTrajectory).x(ti), [], 1);
                yCell{iTrajectory} = reshape(trajectories(iTrajectory).y(ti), [], 1);
                xDotCell{iTrajectory} = reshape(trajectories(iTrajectory).u(ti), [], 1);
                yDotCell{iTrajectory} = reshape(trajectories(iTrajectory).v(ti), [], 1);
            end

            sampleData = struct( ...
                "tCell", {tCell}, ...
                "xCell", {xCell}, ...
                "yCell", {yCell}, ...
                "xDotCell", {xDotCell}, ...
                "yDotCell", {yDotCell}, ...
                "allT", vertcat(tCell{:}), ...
                "allX", vertcat(xCell{:}), ...
                "allY", vertcat(yCell{:}), ...
                "allXDot", vertcat(xDotCell{:}), ...
                "allYDot", vertcat(yDotCell{:}), ...
                "sampleSupportTimes", unique(vertcat(tCell{:}), "sorted"));
        end

        function [fitTrajectories, evaluationTrajectories] = trimmedTrajectoryData()
            [~, t, x, y, fitTrajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            startIndices = [1; 3; 2; 5; 4; 2; 6; 3; 1];
            endIndices = [numel(t) - 2; numel(t); numel(t) - 3; numel(t) - 1; numel(t); numel(t) - 4; numel(t) - 1; numel(t) - 2; numel(t)];

            evaluationTrajectories = TrajectorySpline.empty(0, 1);
            for iTrajectory = 1:size(x, 2)
                indices = startIndices(iTrajectory):2:endIndices(iTrajectory);
                evaluationTrajectories(end + 1, 1) = TrajectorySpline.fromData( ...
                    t(indices), x(indices, iTrajectory), y(indices, iTrajectory), S=3); %#ok<AGROW>
            end
        end

        function psiKnotPoints = zeroPsiKnotPoints()
            psiKnotPoints = {[-1; -1; -1; 1; 1; 1], [-1; -1; -1; 1; 1; 1], [0; 1]};
        end

        function gaugeModeMatrix = constantGaugeMode(psiKnotPoints, psiS, basisSize)
            numSpatialCoefficients = prod(basisSize(1:2));
            if numSpatialCoefficients == 1
                gaugeModeMatrix = 1;
                return
            end

            qDomain = [psiKnotPoints{1}(1), psiKnotPoints{1}(end)];
            rDomain = [psiKnotPoints{2}(1), psiKnotPoints{2}(end)];
            spatialSupportCounts = max(2 * basisSize(1:2), basisSize(1:2) + 2);
            qGaugePoints = linspace(qDomain(1), qDomain(2), spatialSupportCounts(1)).';
            rGaugePoints = linspace(rDomain(1), rDomain(2), spatialSupportCounts(2)).';
            [qGaugeGrid, rGaugeGrid] = ndgrid(qGaugePoints, rGaugePoints);
            spatialPointMatrix = [qGaugeGrid(:), rGaugeGrid(:)];
            spatialBasisMatrix = TensorSpline.matrixForPointMatrix(spatialPointMatrix, knotPoints=psiKnotPoints(1:2), S=psiS(1:2));
            constantGaugeVector = spatialBasisMatrix \ ones(size(qGaugeGrid(:)));
            gaugeModeMatrix = reshape(constantGaugeVector, [], 1);
        end

        function [dof, G, Aeq] = expectedMesoscaleDegreesOfFreedom(fit)
            basisSize = fit.streamfunctionSpline.basisSize;
            numSpatialCoefficients = prod(basisSize(1:2));
            if numSpatialCoefficients == 1
                Gspatial = sparse(1, 0);
            else
                gaugeModeMatrix = GriddedStreamfunctionUnitTests.constantGaugeMode( ...
                    fit.psiKnotPoints, fit.psiS, basisSize);
                Gspatial = sparse(null(reshape(gaugeModeMatrix, [], 1).', "r"));
            end

            G = kron(speye(basisSize(3)), Gspatial);
            Aeq = GriddedStreamfunctionUnitTests.constraintMatrixForFit( ...
                fit.psiKnotPoints, fit.psiS, G, fit.mesoscaleConstraint);
            dof = size(G, 2) - size(Aeq, 1);
        end

        function Aeq = constraintMatrixForFit(psiKnotPoints, psiS, G, mesoscaleConstraint)
            numReducedCoefficients = size(G, 2);
            if mesoscaleConstraint == "none" || numReducedCoefficients == 0
                Aeq = sparse(0, numReducedCoefficients);
                return
            end

            pointMatrix = GriddedStreamfunctionUnitTests.constraintPointMatrix(psiKnotPoints, psiS);
            Bqq = TensorSpline.matrixForPointMatrix(pointMatrix, knotPoints=psiKnotPoints, S=psiS, D=[2 0 0]);
            Brr = TensorSpline.matrixForPointMatrix(pointMatrix, knotPoints=psiKnotPoints, S=psiS, D=[0 2 0]);

            switch mesoscaleConstraint
                case "zeroVorticity"
                    AeqFull = Bqq + Brr;
                case "zeroStrain"
                    Bqr = TensorSpline.matrixForPointMatrix(pointMatrix, knotPoints=psiKnotPoints, S=psiS, D=[1 1 0]);
                    AeqFull = [Bqr; Bqq - Brr];
                otherwise
                    error("GriddedStreamfunctionUnitTests:UnknownMesoscaleConstraint", ...
                        "Unknown mesoscale constraint '%s'.", mesoscaleConstraint);
            end

            Aeq = GriddedStreamfunctionUnitTests.independentRows(AeqFull * G);
        end

        function pointMatrix = constraintPointMatrix(psiKnotPoints, psiS)
            qSupport = GriddedStreamfunctionUnitTests.supportPointsForDimension(psiKnotPoints{1}, psiS(1));
            rSupport = GriddedStreamfunctionUnitTests.supportPointsForDimension(psiKnotPoints{2}, psiS(2));
            tSupport = GriddedStreamfunctionUnitTests.supportPointsForDimension(psiKnotPoints{3}, psiS(3));
            [qGrid, rGrid, tGrid] = ndgrid(qSupport, rSupport, tSupport);
            pointMatrix = [qGrid(:), rGrid(:), tGrid(:)];
        end

        function supportPoints = supportPointsForDimension(knotVector, splineDegree)
            uniqueKnots = unique(knotVector, "sorted");
            intervalStarts = uniqueKnots(1:(end - 1));
            intervalWidths = diff(uniqueKnots);
            numInteriorPoints = splineDegree + 1;
            interiorFractions = (1:numInteriorPoints).' / (numInteriorPoints + 1);
            supportPoints = zeros(numInteriorPoints * nnz(intervalWidths > 0), 1);

            nextIndex = 1;
            for iInterval = 1:numel(intervalStarts)
                if intervalWidths(iInterval) <= 0
                    continue
                end

                intervalPoints = intervalStarts(iInterval) + interiorFractions * intervalWidths(iInterval);
                supportPoints(nextIndex:(nextIndex + numInteriorPoints - 1)) = intervalPoints;
                nextIndex = nextIndex + numInteriorPoints;
            end
        end

        function Aeq = independentRows(Aeq)
            numReducedCoefficients = size(Aeq, 2);
            if isempty(Aeq) || numReducedCoefficients == 0
                Aeq = sparse(0, numReducedCoefficients);
                return
            end

            [~, R, permutation] = qr(full(Aeq.'), "vector");
            diagonal = abs(diag(R));
            if isempty(diagonal)
                rankAeq = 0;
            else
                tolerance = max(size(R)) * eps(class(R)) * max(diagonal);
                rankAeq = nnz(diagonal > tolerance);
            end

            if rankAeq == 0
                Aeq = sparse(0, numReducedCoefficients);
            else
                rowIndices = sort(permutation(1:rankAeq));
                Aeq = sparse(Aeq(rowIndices, :));
            end
        end

        function verifyFiniteReconstruction(testCase, fit)
            testCase.verifyTrue(all(isfinite(fit.backgroundTrajectory.x(fit.fitSupportTimes))))
            testCase.verifyTrue(all(isfinite(fit.backgroundTrajectory.y(fit.fitSupportTimes))))

            for iTrajectory = 1:numel(fit.observedTrajectories)
                trajectory = fit.observedTrajectories(iTrajectory);
                ti = trajectory.t;
                observedX = trajectory.x(ti);
                observedY = trajectory.y(ti);

                background = fit.decomposition.fixedFrame.background(iTrajectory);
                mesoscale = fit.decomposition.fixedFrame.mesoscale(iTrajectory);
                submesoscale = fit.decomposition.fixedFrame.submesoscale(iTrajectory);

                reconstructedX = background.x(ti) + mesoscale.x(ti) + submesoscale.x(ti);
                reconstructedY = background.y(ti) + mesoscale.y(ti) + submesoscale.y(ti);
                componentValues = [ ...
                    background.x(ti); background.y(ti); ...
                    mesoscale.x(ti); mesoscale.y(ti); ...
                    submesoscale.x(ti); submesoscale.y(ti)];

                testCase.verifyTrue(all(isfinite(componentValues)))
                testCase.verifyEqual(reconstructedX, observedX, "AbsTol", 1e-2)
                testCase.verifyEqual(reconstructedY, observedY, "AbsTol", 1e-2)
            end
        end

        function verifyTrajectorySplineEqualAtSamples(testCase, actual, expected, ti)
            testCase.verifyEqual(actual.x(ti), expected.x(ti), "AbsTol", 1e-12)
            testCase.verifyEqual(actual.y(ti), expected.y(ti), "AbsTol", 1e-12)
            testCase.verifyEqual(actual.u(ti), expected.u(ti), "AbsTol", 1e-12)
            testCase.verifyEqual(actual.v(ti), expected.v(ti), "AbsTol", 1e-12)
        end

        function deleteIfExists(path)
            if isfile(path)
                delete(path);
            end
        end
    end
end
