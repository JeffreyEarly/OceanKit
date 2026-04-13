classdef KernelDensityEstimateUnitTests < matlab.unittest.TestCase

    methods (Test)
        function fitOneDimensionalDataMatchesLegacyBandwidth(testCase)
            data = KernelDensityEstimateUnitTests.oneDimensionalData();
            model = KernelDensityEstimate.fromData(data);

            testCase.verifyClass(model, "KernelDensityEstimate")
            testCase.verifyEqual(model.numDimensions, 1)
            testCase.verifyEqual(model.bandwidth, 2.21333986156179, "AbsTol", 1e-12)
        end

        function fitOneDimensionalDataMatchesFrozenBandwidthOnAdditionalDataset(testCase)
            datasets = KernelDensityEstimateUnitTests.oneDimensionalDatasets();
            model = KernelDensityEstimate.fromData(datasets(2).data);

            testCase.verifyEqual(model.bandwidth, datasets(2).bandwidth, "AbsTol", 1e-12)
        end

        function fitOneDimensionalDataIsPermutationInvariant(testCase)
            datasets = KernelDensityEstimateUnitTests.oneDimensionalDatasets();

            for iDataset = 1:numel(datasets)
                data = datasets(iDataset).data;
                baseline = KernelDensityEstimate.fromData(data);
                permuted = KernelDensityEstimate.fromData(data(end:-1:1));

                testCase.verifyEqual(permuted.bandwidth, baseline.bandwidth, "AbsTol", 1e-12)
            end
        end

        function fitOneDimensionalDataIsTranslationInvariant(testCase)
            datasets = KernelDensityEstimateUnitTests.oneDimensionalDatasets();
            translation = 8;

            for iDataset = 1:numel(datasets)
                baseline = KernelDensityEstimate.fromData(datasets(iDataset).data);
                translated = KernelDensityEstimate.fromData(datasets(iDataset).data + translation);

                testCase.verifyEqual(translated.bandwidth, baseline.bandwidth, "AbsTol", 1e-12)
            end
        end

        function fitOneDimensionalDataIsPositiveScaleEquivariant(testCase)
            datasets = KernelDensityEstimateUnitTests.oneDimensionalDatasets();
            scaleFactor = 4;

            for iDataset = 1:numel(datasets)
                baseline = KernelDensityEstimate.fromData(datasets(iDataset).data);
                scaled = KernelDensityEstimate.fromData(datasets(iDataset).data * scaleFactor);

                testCase.verifyEqual(scaled.bandwidth, baseline.bandwidth * scaleFactor, "AbsTol", 1e-10, "RelTol", 1e-10)
            end
        end

        function fitTwoDimensionalDataPreservesLegacyKde2dBandwidthSmokeTest(testCase)
            data = KernelDensityEstimateUnitTests.twoDimensionalData();
            model = KernelDensityEstimate.fromData(data);

            testCase.verifyEqual(model.numDimensions, 2)
            testCase.verifyEqual(model.bandwidth, [1.78508660817316 1.09168358178069], "AbsTol", 1e-12)
        end

        function fitTwoDimensionalDataMatchesReferenceBandwidthOnFrozenDatasets(testCase)
            datasets = KernelDensityEstimateUnitTests.twoDimensionalDatasets();

            for iDataset = 1:numel(datasets)
                model = KernelDensityEstimate.fromData(datasets(iDataset).data);
                diagnostics = referenceBandwidthForTwoDimensionalData(datasets(iDataset).data);

                testCase.verifyEqual(model.bandwidth, diagnostics.bandwidth, "AbsTol", 1e-10, "RelTol", 1e-10)
            end
        end

        function referenceSolverResidualIsSmallOnFrozenDatasets(testCase)
            datasets = KernelDensityEstimateUnitTests.twoDimensionalDatasets();

            for iDataset = 1:numel(datasets)
                diagnostics = referenceBandwidthForTwoDimensionalData(datasets(iDataset).data);

                testCase.verifyLessThanOrEqual(abs(diagnostics.objectiveResidual), 1e-10)
            end
        end

        function fitTwoDimensionalDataIsPermutationInvariant(testCase)
            datasets = KernelDensityEstimateUnitTests.twoDimensionalDatasets();

            for iDataset = 1:numel(datasets)
                data = datasets(iDataset).data;
                baseline = KernelDensityEstimate.fromData(data);
                permuted = KernelDensityEstimate.fromData(data(end:-1:1, :));

                testCase.verifyEqual(permuted.bandwidth, baseline.bandwidth, "AbsTol", 1e-12)
            end
        end

        function fitTwoDimensionalDataIsTranslationInvariant(testCase)
            datasets = KernelDensityEstimateUnitTests.twoDimensionalDatasets();
            translation = [8 -4];

            for iDataset = 1:numel(datasets)
                baseline = KernelDensityEstimate.fromData(datasets(iDataset).data);
                translated = KernelDensityEstimate.fromData(datasets(iDataset).data + translation);

                testCase.verifyEqual(translated.bandwidth, baseline.bandwidth, "AbsTol", 1e-12)
            end
        end

        function fitTwoDimensionalDataIsComponentwiseScaleInvariant(testCase)
            datasets = KernelDensityEstimateUnitTests.twoDimensionalDatasets();
            scaleFactors = [4 0.25];

            for iDataset = 1:numel(datasets)
                baseline = KernelDensityEstimate.fromData(datasets(iDataset).data);
                scaled = KernelDensityEstimate.fromData(datasets(iDataset).data .* scaleFactors);

                testCase.verifyEqual(scaled.bandwidth, baseline.bandwidth .* scaleFactors, "AbsTol", 1e-10, "RelTol", 1e-10)
            end
        end

        function densityAtMatchesDirectGaussianReferenceInOneDimension(testCase)
            data = KernelDensityEstimateUnitTests.oneDimensionalData();
            model = KernelDensityEstimate.fromData(data);
            queryPoints = linspace(-4, 4, 11).';

            expected = KernelDensityEstimateUnitTests.directGaussianDensity1D(model.data, model.bandwidth, queryPoints);

            testCase.verifyEqual(model.densityAt(queryPoints), expected, "AbsTol", 1e-12)
            testCase.verifyEqual(model.densityAt(queryPoints.'), expected.', "AbsTol", 1e-12)
        end

        function densityAtMatchesDirectGaussianReferenceInTwoDimensions(testCase)
            data = KernelDensityEstimateUnitTests.twoDimensionalData();
            model = KernelDensityEstimate.fromData(data);
            queryPoints = [-1.5 -0.8; -0.25 0.5; 1.6 1.2];

            expected = KernelDensityEstimateUnitTests.directGaussianDensity2D(model.data, model.bandwidth, queryPoints);

            testCase.verifyEqual(model.densityAt(queryPoints), expected, "AbsTol", 1e-12)
        end

        function densityOnGridMatchesPointEvaluationAndSupportsArbitraryGridSizes(testCase)
            model1D = KernelDensityEstimate.fromData(KernelDensityEstimateUnitTests.oneDimensionalData());
            [density1D, gridVectors1D] = model1D.densityOnGrid(gridSize=300);

            testCase.verifySize(density1D, [300 1])
            testCase.verifyEqual(gridVectors1D{1}(1), model1D.minimum, "AbsTol", 1e-12)
            testCase.verifyEqual(gridVectors1D{1}(end), model1D.maximum, "AbsTol", 1e-12)
            testCase.verifyEqual(density1D, model1D.densityAt(gridVectors1D{1}), "AbsTol", 1e-12)

            model2D = KernelDensityEstimate.fromData(KernelDensityEstimateUnitTests.twoDimensionalData());
            [density2D, gridVectors2D] = model2D.densityOnGrid(gridSize=[192 320]);
            [X, Y] = ndgrid(gridVectors2D{1}, gridVectors2D{2});

            testCase.verifySize(density2D, [192 320])
            testCase.verifyEqual(density2D, reshape(model2D.densityAt([X(:), Y(:)]), size(X)), "AbsTol", 1e-12)
        end

        function cdfOnGridMatchesDirectGaussianReference(testCase)
            data = KernelDensityEstimateUnitTests.oneDimensionalData();
            model = KernelDensityEstimate.fromData(data);
            [cdfValues, gridVectors] = model.cdfOnGrid(gridSize=301);

            expected = KernelDensityEstimateUnitTests.directGaussianCdf1D(model.data, model.bandwidth, gridVectors{1});

            testCase.verifyEqual(cdfValues, expected, "AbsTol", 1e-12)
            testCase.verifyGreaterThanOrEqual(min(cdfValues), 0)
            testCase.verifyLessThanOrEqual(max(cdfValues), 1)
        end

        function bandwidthGridSizeRoundsUpInternally(testCase)
            data = KernelDensityEstimateUnitTests.oneDimensionalData();
            roundedModel = KernelDensityEstimate.fromData(data, bandwidthGridSize=300);
            exactModel = KernelDensityEstimate.fromData(data, bandwidthGridSize=512);

            testCase.verifyEqual(roundedModel.bandwidth, exactModel.bandwidth, "AbsTol", 1e-12)
        end

        function twoDimensionalBandwidthGridSizeRoundsUpInternally(testCase)
            datasets = KernelDensityEstimateUnitTests.twoDimensionalDatasets();

            for iDataset = 1:numel(datasets)
                roundedModel = KernelDensityEstimate.fromData(datasets(iDataset).data, bandwidthGridSize=300);
                exactModel = KernelDensityEstimate.fromData(datasets(iDataset).data, bandwidthGridSize=512);

                testCase.verifyEqual(roundedModel.bandwidth, exactModel.bandwidth, "AbsTol", 1e-12)
            end
        end

        function invalidDimensionsAndBoundsRaiseStructuredErrors(testCase)
            data = KernelDensityEstimateUnitTests.oneDimensionalData();

            testCase.verifyError(@() KernelDensityEstimate.fromData(rand(5, 3)), "KernelDensityEstimate:UnsupportedDimension")
            testCase.verifyError(@() KernelDensityEstimate.fromData(data, minimum=1, maximum=0), "KernelDensityEstimate:InvalidBounds")
            testCase.verifyError(@() KernelDensityEstimate.fromData(data, minimum=[-1 -1], maximum=[1 1]), "KernelDensityEstimate:InvalidBounds")
        end

        function densityLevelForCdfReturnsTargetMassThresholds(testCase)
            data = KernelDensityEstimateUnitTests.twoDimensionalData();
            model = KernelDensityEstimate.fromData(data);
            [density, gridVectors] = model.densityOnGrid(gridSize=[64 64]);
            pctTarget = [0.1; 0.5; 0.8];

            levels = DensityLevelForCDF(gridVectors, density, pctTarget);
            enclosedMass = KernelDensityEstimateUnitTests.enclosedMassForDensityLevel(gridVectors, density, levels);

            testCase.verifySize(levels, size(pctTarget))
            testCase.verifyTrue(all(isfinite(levels)))
            testCase.verifyTrue(all(diff(levels) <= 0))
            testCase.verifyLessThanOrEqual(max(abs(enclosedMass - pctTarget)), 0.08)
        end

        function densityLevelForCdfHandlesDuplicateEnclosedMassSamples(testCase)
            xVector = reshape(linspace(-5, 5, 128), [], 1);
            yVector = reshape(linspace(-4, 4, 128), [], 1);
            [X, Y] = ndgrid(xVector, yVector);
            density = exp(-((X - 1.5).^2 + Y.^2)) + exp(-((X + 1.5).^2 + Y.^2));
            density = density / trapz(yVector, trapz(xVector, density, 1), 2);
            pctTarget = [0.05; 0.1; 0.15];

            actual = DensityLevelForCDF({xVector, yVector}, density, pctTarget);

            testCase.verifySize(actual, size(pctTarget))
            testCase.verifyTrue(all(isfinite(actual)))
            testCase.verifyTrue(all(diff(actual) <= 0))
        end

        function planarStatisticsFromDataBuildsPlanarSummary(testCase)
            data = KernelDensityEstimateUnitTests.twoDimensionalData();
            statistics = KernelDensityEstimate.planarStatisticsFromData(data, referencePoint=[0.5 -0.25], gridSize=[96 112], summaryMass=0.5);

            testCase.verifyEqual(size(statistics.modePoint), [1 2])
            testCase.verifyEqual(size(statistics.referencePoint), [1 2])
            testCase.verifyEqual(size(statistics.density), [96 112])
            testCase.verifyEqual(statistics.referencePoint, [0.5 -0.25], "AbsTol", 1e-12)
            testCase.verifyTrue(all(isfinite(statistics.modePoint)))
            testCase.verifyTrue(all(isfinite(statistics.radiusBounds)))
            testCase.verifyTrue(all(isfinite(statistics.angleBounds)))
            testCase.verifyTrue(all(statistics.bounds.maximum > statistics.bounds.minimum))
            testCase.verifyNotEmpty(statistics.summaryContourVertices)
            testCase.verifyTrue(inpolygon( ...
                statistics.modePoint(1), ...
                statistics.modePoint(2), ...
                statistics.summaryContourVertices(:, 1), ...
                statistics.summaryContourVertices(:, 2)))
        end

        function planarStatisticsFromDataSelectsDisconnectedSummaryContour(testCase)
            datasets = KernelDensityEstimateUnitTests.twoDimensionalDatasets();
            statistics = KernelDensityEstimate.planarStatisticsFromData(datasets(2).data, gridSize=[128 128], summaryMass=0.15);

            testCase.verifyNotEmpty(statistics.summaryContourVertices)
            testCase.verifyTrue(inpolygon( ...
                statistics.modePoint(1), ...
                statistics.modePoint(2), ...
                statistics.summaryContourVertices(:, 1), ...
                statistics.summaryContourVertices(:, 2)))
            testCase.verifyLessThan(range(statistics.summaryContourVertices(:, 1)), 4)
        end

        function planarStatisticsFromDataHandlesMissingSummaryContour(testCase)
            data = KernelDensityEstimateUnitTests.twoDimensionalData();
            statistics = KernelDensityEstimate.planarStatisticsFromData(data, summaryMass=0);

            testCase.verifyTrue(isnan(statistics.summaryLevel))
            testCase.verifyEmpty(statistics.summaryContourVertices)
            testCase.verifyTrue(all(isnan(statistics.radiusBounds)))
            testCase.verifyTrue(all(isnan(statistics.angleBounds)))
            testCase.verifyTrue(all(isfinite(statistics.modePoint)))
        end

        function planarStatisticsFromDataValidatesPlanarOptionShapes(testCase)
            data = KernelDensityEstimateUnitTests.twoDimensionalData();

            testCase.verifyError( ...
                @() KernelDensityEstimate.planarStatisticsFromData(data, referencePoint=[1 2 3]), ...
                "KernelDensityEstimate:InvalidPoint")
            testCase.verifyError( ...
                @() KernelDensityEstimate.planarStatisticsFromData(data, minimum=[-1 -1 -1], maximum=[1 1 1]), ...
                "KernelDensityEstimate:InvalidBounds")
        end

        function polarSummaryFromPlanarStatisticsSupportsAngleReference(testCase)
            data = KernelDensityEstimateUnitTests.twoDimensionalData();
            statistics = KernelDensityEstimate.planarStatisticsFromData(data, referencePoint=[0.5 -0.25], summaryMass=0.5);

            defaultPolar = KernelDensityEstimate.polarSummaryFromPlanarStatistics(statistics);
            shiftedPolar = KernelDensityEstimate.polarSummaryFromPlanarStatistics(statistics, angleReference=defaultPolar.modeAngle + 2*pi);

            testCase.verifyClass(defaultPolar.containsOrigin, "logical")
            testCase.verifySize(defaultPolar.containsOrigin, [1 1])
            testCase.verifyEqual(shiftedPolar.modeRadius, defaultPolar.modeRadius, "AbsTol", 1e-12)
            testCase.verifyEqual(shiftedPolar.modeAngle, defaultPolar.modeAngle + 2*pi, "AbsTol", 1e-12)
            testCase.verifyEqual(shiftedPolar.radiusBounds, defaultPolar.radiusBounds, "AbsTol", 1e-12)
            if all(isfinite(defaultPolar.angleBounds))
                testCase.verifyEqual(shiftedPolar.angleBounds, defaultPolar.angleBounds + 2*pi, "AbsTol", 1e-12)
            end
        end

        function polarSummaryFromPlanarStatisticsMarksOriginInclusiveContours(testCase)
            statistics = KernelDensityEstimateUnitTests.manualPlanarStatistics( ...
                [1 1], ...
                KernelDensityEstimateUnitTests.circleVertices([0 0], 2, 48), ...
                referencePoint=[0.5 -0.25], ...
                radiusBounds=[0.5 2.0], ...
                angleBounds=[0 1]);

            polarSummary = KernelDensityEstimate.polarSummaryFromPlanarStatistics(statistics);

            testCase.verifyTrue(polarSummary.containsOrigin)
            testCase.verifyEqual(polarSummary.modeRadius, sqrt(2), "AbsTol", 1e-12)
            testCase.verifyEqual(polarSummary.radiusBounds, [0 2], "AbsTol", 1e-12)
            testCase.verifyEqual(polarSummary.radiusErrors, [sqrt(2), 2 - sqrt(2)], "AbsTol", 1e-12)
            testCase.verifyTrue(all(isnan(polarSummary.angleBounds)))
            testCase.verifyTrue(all(isnan(polarSummary.angleErrors)))
        end

        function strainSummaryFromPlanarStatisticsRecoversSigmaThetaOnTightCloud(testCase)
            sigma = 4;
            theta = pi/6;
            center = [sigma * cos(2 * theta), sigma * sin(2 * theta)];
            data = center + [ ...
                0.00 0.00; ...
                0.08 -0.04; ...
                -0.06 0.05; ...
                0.04 0.03; ...
                -0.05 -0.06; ...
                0.03 0.02];
            statistics = KernelDensityEstimate.planarStatisticsFromData( ...
                data, ...
                minimum=center - [0.75 0.75], ...
                maximum=center + [0.75 0.75], ...
                gridSize=[128 128], ...
                summaryMass=0.5);

            strainSummary = KernelDensityEstimate.strainSummaryFromPlanarStatistics(statistics);

            testCase.verifyFalse(strainSummary.containsZero)
            testCase.verifyLessThanOrEqual(abs(strainSummary.modeSigma - sigma), 0.2)
            testCase.verifyLessThanOrEqual(abs(strainSummary.modeTheta - theta), 0.08)
            testCase.verifyTrue(all(isfinite(strainSummary.sigmaBounds)))
            testCase.verifyTrue(all(isfinite(strainSummary.thetaBounds)))
        end

        function strainSummaryFromPlanarStatisticsSupportsThetaReference(testCase)
            statistics = KernelDensityEstimateUnitTests.manualPlanarStatistics( ...
                [2 * cos(pi/3), 2 * sin(pi/3)], ...
                KernelDensityEstimateUnitTests.circleVertices([1 0], 0.5, 48), ...
                radiusBounds=[1.5 2.5], ...
                angleBounds=[pi/6 pi/2]);

            defaultSummary = KernelDensityEstimate.strainSummaryFromPlanarStatistics(statistics);
            shiftedSummary = KernelDensityEstimate.strainSummaryFromPlanarStatistics( ...
                statistics, ...
                thetaReference=defaultSummary.modeTheta + pi);

            testCase.verifyEqual(shiftedSummary.modeSigma, defaultSummary.modeSigma, "AbsTol", 1e-12)
            testCase.verifyEqual(shiftedSummary.modeTheta, defaultSummary.modeTheta + pi, "AbsTol", 1e-12)
            testCase.verifyEqual(shiftedSummary.sigmaBounds, defaultSummary.sigmaBounds, "AbsTol", 1e-12)
            testCase.verifyEqual(shiftedSummary.thetaBounds, defaultSummary.thetaBounds + pi, "AbsTol", 1e-12)
        end

        function strainSummaryFromPlanarStatisticsKeepsOppositeModesDistinct(testCase)
            statisticsA = KernelDensityEstimateUnitTests.manualPlanarStatistics( ...
                [2 0], ...
                KernelDensityEstimateUnitTests.circleVertices([2 0], 0.4, 48), ...
                radiusBounds=[1.6 2.4], ...
                angleBounds=[-0.2 0.2]);
            statisticsB = KernelDensityEstimateUnitTests.manualPlanarStatistics( ...
                [-2 0], ...
                KernelDensityEstimateUnitTests.circleVertices([-2 0], 0.4, 48), ...
                radiusBounds=[1.6 2.4], ...
                angleBounds=[pi - 0.2 pi + 0.2]);

            summaryA = KernelDensityEstimate.strainSummaryFromPlanarStatistics(statisticsA);
            summaryB = KernelDensityEstimate.strainSummaryFromPlanarStatistics(statisticsB);

            testCase.verifyEqual(summaryB.modeTheta - summaryA.modeTheta, pi/2, "AbsTol", 1e-12)
        end

        function strainSummaryFromPlanarStatisticsPreservesModeWhenContourContainsOrigin(testCase)
            statistics = KernelDensityEstimateUnitTests.manualPlanarStatistics( ...
                [1 1], ...
                KernelDensityEstimateUnitTests.circleVertices([0 0], 2, 48), ...
                radiusBounds=[0.5 2.0], ...
                angleBounds=[0 1]);

            strainSummary = KernelDensityEstimate.strainSummaryFromPlanarStatistics(statistics);

            testCase.verifyTrue(strainSummary.containsZero)
            testCase.verifyEqual(strainSummary.modeSigma, sqrt(2), "AbsTol", 1e-12)
            testCase.verifyEqual(strainSummary.sigmaBounds, [0 2], "AbsTol", 1e-12)
            testCase.verifyEqual(strainSummary.modeTheta, pi/8, "AbsTol", 1e-12)
            testCase.verifyTrue(all(isnan(strainSummary.thetaBounds)))
            testCase.verifyTrue(all(isnan(strainSummary.thetaErrors)))
        end

        function planarHelpersValidateStatisticsShapesAtThePublicBoundary(testCase)
            data = KernelDensityEstimateUnitTests.twoDimensionalData();
            statistics = KernelDensityEstimate.planarStatisticsFromData(data, summaryMass=0.5);
            invalidPlotStatistics = statistics;
            invalidPolarStatistics = statistics;
            invalidPlotStatistics.angleBounds = [statistics.angleBounds, statistics.angleBounds(end)];
            invalidPolarStatistics.radiusBounds = [statistics.radiusBounds, statistics.radiusBounds(end)];

            figureHandle = figure(Visible="off");
            testCase.addTeardown(@close, figureHandle);
            ax = axes(figureHandle);

            testCase.verifyError(@() KernelDensityEstimate.plotPlanarStatistics(ax, invalidPlotStatistics), ...
                "KernelDensityEstimate:InvalidPlanarStatistics")
            testCase.verifyError(@() KernelDensityEstimate.polarSummaryFromPlanarStatistics(invalidPolarStatistics), ...
                "KernelDensityEstimate:InvalidPlanarStatistics")
        end

        function plotPlanarStatisticsDrawsPlanarSummary(testCase)
            data = KernelDensityEstimateUnitTests.twoDimensionalData();
            statistics = KernelDensityEstimate.planarStatisticsFromData(data, referencePoint=[0.5 -0.25], gridSize=[64 64], summaryMass=0.5);

            figureHandle = figure(Visible="off");
            testCase.addTeardown(@close, figureHandle);
            ax = axes(figureHandle);

            KernelDensityEstimate.plotPlanarStatistics(ax, statistics, ...
                samples=data, ...
                contourMasses=[0.2; 0.5; 0.8], ...
                showWedge=true, ...
                ringRadii=[0.5; 1.0]);

            testCase.verifyEqual(ax.XLim, [statistics.bounds.minimum(1), statistics.bounds.maximum(1)], "AbsTol", 1e-12)
            testCase.verifyEqual(ax.YLim, [statistics.bounds.minimum(2), statistics.bounds.maximum(2)], "AbsTol", 1e-12)
            testCase.verifyGreaterThan(numel(ax.Children), 0)
        end

        function plotPlanarStatisticsSkipsWedgeWhenAngleBoundsAreUndefined(testCase)
            data = KernelDensityEstimateUnitTests.twoDimensionalData();
            statistics = KernelDensityEstimate.planarStatisticsFromData(data, referencePoint=[0.5 -0.25], gridSize=[64 64], summaryMass=0.5);
            statistics.angleBounds = [NaN NaN];

            figureHandle = figure(Visible="off");
            testCase.addTeardown(@close, figureHandle);
            ax = axes(figureHandle);

            KernelDensityEstimate.plotPlanarStatistics(ax, statistics, ...
                samples=data, ...
                contourMasses=[0.2; 0.5; 0.8], ...
                showWedge=true);

            testCase.verifyGreaterThan(numel(ax.Children), 0)
        end
    end

    methods (Static, Access = private)
        function data = oneDimensionalData()
            datasets = KernelDensityEstimateUnitTests.oneDimensionalDatasets();
            data = datasets(1).data;
        end

        function datasets = oneDimensionalDatasets()
            datasets = [ ...
                struct( ...
                    "name", "Original sparse sample", ...
                    "data", [-3.1; -1.2; -0.8; 0.1; 0.4; 0.9; 1.7; 2.4; 3.0], ...
                    "bandwidth", 2.21333986156179), ...
                struct( ...
                    "name", "Bimodal cloud with repeated sample", ...
                    "data", [-5.0; -4.5; -4.0; -3.0; 2.5; 3.0; 3.0; 4.0; 4.5; 5.0], ...
                    "bandwidth", 1.73785891998607)];
        end

        function data = twoDimensionalData()
            data = [-2.0 -1.1; -1.4 -0.2; -0.8 0.4; -0.2 0.9; 0.5 -0.6; 1.0 0.2; 1.7 1.1; 2.4 1.8];
        end

        function datasets = twoDimensionalDatasets()
            datasets = [ ...
                struct( ...
                    "name", "Compact anisotropic cloud", ...
                    "data", [-2.00 -1.00; -1.50 -0.25; -1.00 0.50; -0.25 1.00; 0.50 -0.50; 1.00 0.25; 1.75 1.00; 2.50 1.75]), ...
                struct( ...
                    "name", "Well-separated axis-aligned bimodal cloud", ...
                    "data", [-5.0 -1.50; -4.5 -1.00; -4.0 -1.25; -3.5 -0.50; -3.0 -1.00; 2.5 1.50; 3.0 2.00; 3.5 1.75; 4.0 2.25; 4.5 2.00]), ...
                struct( ...
                    "name", "Strongly anisotropic single cloud", ...
                    "data", [-6.00 0.12500; -4.75 -0.12500; -3.50 0.06250; -2.25 -0.06250; -1.00 0.03125; 0.25 -0.03125; 1.50 0.09375; 2.75 -0.09375; 4.00 0.046875; 5.50 -0.046875])];
        end

        function density = directGaussianDensity1D(data, bandwidth, queryPoints)
            queryPoints = reshape(queryPoints, [], 1);
            density = mean(exp(-0.5 * ((queryPoints - data.')/bandwidth).^2), 2)/(bandwidth * sqrt(2*pi));
        end

        function density = directGaussianDensity2D(data, bandwidth, queryPoints)
            density = mean(exp(-0.5 * (((queryPoints(:, 1) - data(:, 1).')/bandwidth(1)).^2 + ((queryPoints(:, 2) - data(:, 2).')/bandwidth(2)).^2)), 2) ...
                /(2*pi * bandwidth(1) * bandwidth(2));
        end

        function cdfValues = directGaussianCdf1D(data, bandwidth, queryPoints)
            queryPoints = reshape(queryPoints, [], 1);
            cdfValues = mean(0.5 * erfc(-(queryPoints - data.')/(bandwidth * sqrt(2))), 2);
        end

        function enclosedMass = enclosedMassForDensityLevel(gridVectors, density, levels)
            xVector = reshape(gridVectors{1}, [], 1);
            yVector = reshape(gridVectors{2}, [], 1);
            enclosedMass = NaN(size(levels));

            for iLevel = 1:numel(levels)
                if ~isfinite(levels(iLevel))
                    continue
                end

                mask = density >= levels(iLevel);
                enclosedMass(iLevel) = trapz(yVector, trapz(xVector, density .* mask, 1), 2);
            end
        end

        function statistics = manualPlanarStatistics(modePoint, summaryContourVertices, options)
            arguments
                modePoint (1,2) double {mustBeReal, mustBeFinite}
                summaryContourVertices (:,2) double {mustBeReal, mustBeFinite}
                options.referencePoint (1,2) double {mustBeReal} = [NaN NaN]
                options.originPoint (1,2) double {mustBeReal, mustBeFinite} = [0 0]
                options.radiusBounds (1,2) double {mustBeReal}
                options.angleBounds (1,2) double {mustBeReal}
            end

            modeOffset = modePoint - options.originPoint;
            statistics = struct( ...
                "modeRadius", hypot(modeOffset(1), modeOffset(2)), ...
                "modeAngle", atan2(modeOffset(2), modeOffset(1)), ...
                "radiusBounds", options.radiusBounds, ...
                "angleBounds", options.angleBounds, ...
                "referencePoint", options.referencePoint, ...
                "originPoint", options.originPoint, ...
                "summaryContourVertices", summaryContourVertices);
        end

        function vertices = circleVertices(center, radius, nVertices)
            angles = linspace(0, 2*pi, nVertices + 1).';
            vertices = center + radius * [cos(angles), sin(angles)];
        end
    end
end
