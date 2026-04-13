classdef GriddedStreamfunctionBootstrapUnitTests < matlab.unittest.TestCase

    methods (Test)
        function constructorStoresOutputsAndReproducesSeed(testCase)
            [~, t, ~, ~, trajectories] = GriddedStreamfunctionBootstrapUnitTests.synchronousLinearFieldData();

            bootstrapA = GriddedStreamfunctionBootstrap.fromTrajectories(trajectories, nBootstraps=6, randomSeed=17);
            bootstrapB = GriddedStreamfunctionBootstrap.fromTrajectories(trajectories, nBootstraps=6, randomSeed=17);

            testCase.verifyClass(bootstrapA, "GriddedStreamfunctionBootstrap")
            testCase.verifyTrue(isa(bootstrapA, "handle"))
            testCase.verifyClass(bootstrapA.fullFit, "GriddedStreamfunction")
            testCase.verifyEqual(bootstrapA.nBootstraps, 6)
            testCase.verifyEqual(bootstrapA.randomSeed, 17)
            testCase.verifyEqual(size(bootstrapA.bootstrapIndices), [6 numel(trajectories)])
            testCase.verifyEqual(bootstrapA.queryTimes, t, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrapA.scoreTimes, t, "AbsTol", 1e-12)
            testCase.verifySize(bootstrapA.summary.uCenter, [numel(t), 6])
            testCase.verifySize(bootstrapA.summary.sigma_n, [numel(t), 6])
            testCase.verifySize(bootstrapA.bootstrapKappa, [1 6])
            testCase.verifySize(bootstrapA.bootstrapCoherence, [1 6])
            testCase.verifyEqual(bootstrapA.mesoscaleDegreesOfFreedom, bootstrapA.fullFit.mesoscaleDegreesOfFreedom)
            testCase.verifyTrue(all(isfinite(bootstrapA.scores.uv), "all"))
            testCase.verifyTrue(all(isfinite(bootstrapA.scores.strain), "all"))
            testCase.verifyTrue(all(isfinite(bootstrapA.scores.zeta), "all"))
            testCase.verifyTrue(all(isfinite(bootstrapA.scores.joint), "all"))
            testCase.verifyGreaterThanOrEqual(min(bootstrapA.bootstrapIndices, [], "all"), 1)
            testCase.verifyLessThanOrEqual(max(bootstrapA.bootstrapIndices, [], "all"), numel(trajectories))

            testCase.verifyEqual(bootstrapA.bootstrapIndices, bootstrapB.bootstrapIndices)
            testCase.verifyEqual(bootstrapA.summary.uCenter, bootstrapB.summary.uCenter, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrapA.summary.vCenter, bootstrapB.summary.vCenter, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrapA.summary.sigma_n, bootstrapB.summary.sigma_n, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrapA.summary.sigma_s, bootstrapB.summary.sigma_s, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrapA.summary.zeta, bootstrapB.summary.zeta, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrapA.bootstrapKappa, bootstrapB.bootstrapKappa, "AbsTol", 1e-12)
            testCase.verifyTrue(isequaln(bootstrapA.bootstrapCoherence, bootstrapB.bootstrapCoherence))
            testCase.verifyEqual(bootstrapA.mesoscaleDegreesOfFreedom, bootstrapB.mesoscaleDegreesOfFreedom)
            testCase.verifyEqual(bootstrapA.fullFitKappa, bootstrapB.fullFitKappa, "AbsTol", 1e-12)
            testCase.verifyTrue(isequaln(bootstrapA.fullFitCoherence, bootstrapB.fullFitCoherence))
            testCase.verifyEqual(bootstrapA.mesoscaleDegreesOfFreedom, bootstrapB.mesoscaleDegreesOfFreedom)
            testCase.verifyEqual(bootstrapA.scores.uv, bootstrapB.scores.uv, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrapA.scores.strain, bootstrapB.scores.strain, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrapA.scores.zeta, bootstrapB.scores.zeta, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrapA.scores.joint, bootstrapB.scores.joint, "AbsTol", 1e-12)

            for iBootstrap = 1:bootstrapA.nBootstraps
                testCase.verifyEqual(bootstrapA.bootstrapMetadata.fastKnotPoints{iBootstrap}, ...
                    bootstrapB.bootstrapMetadata.fastKnotPoints{iBootstrap}, "AbsTol", 1e-12)
                for iDim = 1:3
                    testCase.verifyEqual(bootstrapA.bootstrapMetadata.psiKnotPoints{iBootstrap}{iDim}, ...
                        bootstrapB.bootstrapMetadata.psiKnotPoints{iBootstrap}{iDim}, "AbsTol", 1e-12)
                end
            end
        end

        function persistenceRoundTripPreservesBootstrapState(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionBootstrapUnitTests.synchronousLinearFieldData();
            bootstrap = GriddedStreamfunctionBootstrap.fromTrajectories(trajectories, nBootstraps=4, randomSeed=12);
            path = [tempname, '.nc'];
            cleanup = onCleanup(@() GriddedStreamfunctionBootstrapUnitTests.deleteFileIfExists(path));

            bootstrap.writeToFile(path, shouldOverwriteExisting=true);
            ncfile = NetCDFFile(path, shouldReadOnly=true);
            cleanupFile = onCleanup(@() ncfile.close());
            group = ncfile;
            if ncfile.hasGroupWithName('GriddedStreamfunctionBootstrap')
                group = ncfile.groupWithName('GriddedStreamfunctionBootstrap');
            end
            testCase.verifyFalse(group.hasVariableWithName('bootstrapKappaCache'))
            testCase.verifyFalse(group.hasVariableWithName('bootstrapCoherenceCache'))
            restored = GriddedStreamfunctionBootstrap.fromFile(path);

            testCase.verifyEqual(restored.nBootstraps, bootstrap.nBootstraps)
            testCase.verifyEqual(restored.randomSeed, bootstrap.randomSeed)
            testCase.verifyEqual(restored.queryTimes, bootstrap.queryTimes, "AbsTol", 1e-12)
            testCase.verifyEqual(restored.scoreTimes, bootstrap.scoreTimes, "AbsTol", 1e-12)
            testCase.verifyEqual(restored.bootstrapIndices, bootstrap.bootstrapIndices)
            testCase.verifyEqual(restored.summary.uCenter, bootstrap.summary.uCenter, "AbsTol", 1e-12)
            testCase.verifyEqual(restored.summary.vCenter, bootstrap.summary.vCenter, "AbsTol", 1e-12)
            testCase.verifyEqual(restored.summary.sigma_n, bootstrap.summary.sigma_n, "AbsTol", 1e-12)
            testCase.verifyEqual(restored.summary.sigma_s, bootstrap.summary.sigma_s, "AbsTol", 1e-12)
            testCase.verifyEqual(restored.summary.zeta, bootstrap.summary.zeta, "AbsTol", 1e-12)
            testCase.verifyEqual(restored.fullSummary.uCenter, bootstrap.fullSummary.uCenter, "AbsTol", 1e-12)
            testCase.verifyEqual(restored.fullFitKappa, bootstrap.fullFitKappa, "AbsTol", 1e-12)
            testCase.verifyTrue(isequaln(restored.fullFitCoherence, bootstrap.fullFitCoherence))
            testCase.verifyEqual(restored.mesoscaleDegreesOfFreedom, bootstrap.mesoscaleDegreesOfFreedom)
            testCase.verifyEqual(restored.bestFitKappa, bootstrap.bestFitKappa, "AbsTol", 1e-12)
            testCase.verifyTrue(isequaln(restored.bestFitCoherence, bootstrap.bestFitCoherence))
            testCase.verifyEqual(restored.scores.uv, bootstrap.scores.uv, "AbsTol", 1e-12)
            testCase.verifyEqual(restored.scores.strain, bootstrap.scores.strain, "AbsTol", 1e-12)
            testCase.verifyEqual(restored.scores.zeta, bootstrap.scores.zeta, "AbsTol", 1e-12)
            testCase.verifyEqual(restored.scores.joint, bootstrap.scores.joint, "AbsTol", 1e-12)
            testCase.verifyEqual(restored.bestBootstrapIndex(), bootstrap.bestBootstrapIndex())

            restoredFit = restored.fitForBootstrap(2);
            referenceFit = bootstrap.fitForBootstrap(2);
            [restoredSummary, restoredDiagnostics] = GriddedStreamfunctionBootstrapUnitTests.summaryFromFit(restoredFit, restored.queryTimes);
            [referenceSummary, referenceDiagnostics] = GriddedStreamfunctionBootstrapUnitTests.summaryFromFit(referenceFit, bootstrap.queryTimes);

            testCase.verifyEqual(restoredSummary.uCenter, referenceSummary.uCenter, "AbsTol", 1e-12)
            testCase.verifyEqual(restoredSummary.vCenter, referenceSummary.vCenter, "AbsTol", 1e-12)
            testCase.verifyEqual(restoredSummary.sigma_n, referenceSummary.sigma_n, "AbsTol", 1e-12)
            testCase.verifyEqual(restoredSummary.sigma_s, referenceSummary.sigma_s, "AbsTol", 1e-12)
            testCase.verifyEqual(restoredSummary.zeta, referenceSummary.zeta, "AbsTol", 1e-12)
            testCase.verifyEqual(restoredDiagnostics.kappa, referenceDiagnostics.kappa, "AbsTol", 1e-12)
        end

        function persistenceRoundTripPreservesMaterializedBootstrapDiagnostics(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionBootstrapUnitTests.synchronousLinearFieldData();
            bootstrap = GriddedStreamfunctionBootstrap.fromTrajectories(trajectories, nBootstraps=4, randomSeed=12);
            bootstrapKappa = bootstrap.bootstrapKappa;
            bootstrapCoherence = bootstrap.bootstrapCoherence;
            path = [tempname, '.nc'];
            cleanup = onCleanup(@() GriddedStreamfunctionBootstrapUnitTests.deleteFileIfExists(path));

            bootstrap.writeToFile(path, shouldOverwriteExisting=true);
            ncfile = NetCDFFile(path, shouldReadOnly=true);
            cleanupFile = onCleanup(@() ncfile.close());
            group = ncfile;
            if ncfile.hasGroupWithName('GriddedStreamfunctionBootstrap')
                group = ncfile.groupWithName('GriddedStreamfunctionBootstrap');
            end
            testCase.verifyTrue(group.hasVariableWithName('bootstrapKappaCache'))
            testCase.verifyTrue(group.hasVariableWithName('bootstrapCoherenceCache'))

            restored = GriddedStreamfunctionBootstrap.fromFile(path);
            testCase.verifyEqual(restored.bootstrapKappa, bootstrapKappa, "AbsTol", 1e-12)
            testCase.verifyTrue(isequaln(restored.bootstrapCoherence, bootstrapCoherence))
        end

        function caseStudyBootstrapCacheReloadsPersistentRestart(testCase)
            [~, t, ~, ~, trajectories] = GriddedStreamfunctionBootstrapUnitTests.synchronousLinearFieldData();
            cacheDirectory = tempname;
            mkdir(cacheDirectory);
            cleanup = onCleanup(@() GriddedStreamfunctionBootstrapUnitTests.deleteDirectoryIfExists(cacheDirectory));

            siteNumber = 1;
            nBootstraps = 4;
            randomSeed = 3;
            scoreStride = 2;
            psiS = [2 2 4];
            fastS = 3;
            mesoscaleConstraint = "none";
            queryTimes = t;
            psiSTag = join(string(psiS), "-");
            cacheFilename = "Rho" + siteNumber + ...
                "GriddedStreamfunctionBootstrapFits" + nBootstraps + ...
                "_seed" + randomSeed + ...
                "_stride" + scoreStride + ...
                "_fastS" + fastS + ...
                "_psiS" + psiSTag + ...
                "_mesoscaleConstraint-" + mesoscaleConstraint + ".nc";
            cachePath = fullfile(cacheDirectory, cacheFilename);

            bootstrapA = GriddedStreamfunctionBootstrap.fromTrajectories( ...
                trajectories, ...
                nBootstraps=nBootstraps, ...
                randomSeed=randomSeed, ...
                queryTimes=queryTimes, ...
                scoreStride=scoreStride, ...
                psiS=psiS, ...
                fastS=fastS, ...
                mesoscaleConstraint=mesoscaleConstraint);
            bootstrapA.writeToFile(cachePath, shouldOverwriteExisting=true);
            bootstrapB = GriddedStreamfunctionBootstrap.fromFile(cachePath);

            testCase.verifyTrue(isfile(cachePath))
            testCase.verifyEqual(bootstrapB.queryTimes, queryTimes, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrapB.scoreTimes, queryTimes(1:scoreStride:end), "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrapB.bootstrapIndices, bootstrapA.bootstrapIndices)
            testCase.verifyEqual(bootstrapB.summary.uCenter, bootstrapA.summary.uCenter, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrapB.summary.vCenter, bootstrapA.summary.vCenter, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrapB.summary.sigma_n, bootstrapA.summary.sigma_n, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrapB.summary.sigma_s, bootstrapA.summary.sigma_s, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrapB.summary.zeta, bootstrapA.summary.zeta, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrapB.scores.joint, bootstrapA.scores.joint, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrapB.fullFitKappa, bootstrapA.fullFitKappa, "AbsTol", 1e-12)
            testCase.verifyTrue(isequaln(bootstrapB.fullFitCoherence, bootstrapA.fullFitCoherence))
            testCase.verifyEqual(bootstrapB.mesoscaleDegreesOfFreedom, bootstrapA.mesoscaleDegreesOfFreedom)
        end

        function fullSummaryMatchesSynchronousLinearFieldTruth(testCase)
            [model, t, ~, ~, trajectories] = GriddedStreamfunctionBootstrapUnitTests.synchronousLinearFieldData();

            bootstrap = GriddedStreamfunctionBootstrap.fromTrajectories(trajectories, nBootstraps=4, randomSeed=5);

            testCase.verifyLessThanOrEqual(max(abs(bootstrap.fullSummary.uCenter - model.u0)), 2e-3)
            testCase.verifyLessThanOrEqual(max(abs(bootstrap.fullSummary.vCenter - model.v0)), 5e-3)
            testCase.verifyLessThanOrEqual(max(abs(bootstrap.fullSummary.sigma_n - model.sigma_n)), 1e-12)
            testCase.verifyLessThanOrEqual(max(abs(bootstrap.fullSummary.sigma_s - model.sigma_s)), 1e-12)
            testCase.verifyLessThanOrEqual(max(abs(bootstrap.fullSummary.zeta - model.zeta)), 1e-12)
            testCase.verifyEqual(bootstrap.queryTimes, t, "AbsTol", 1e-12)
        end

        function bestFitAndStoredReplicateReconstructExactly(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionBootstrapUnitTests.synchronousLinearFieldData();

            bootstrap = GriddedStreamfunctionBootstrap.fromTrajectories(trajectories, nBootstraps=6, randomSeed=9);
            iBest = bootstrap.bestBootstrapIndex();
            fit = bootstrap.bestFit();
            [summary, diagnostics] = GriddedStreamfunctionBootstrapUnitTests.summaryFromFit(fit, bootstrap.queryTimes);

            testCase.verifyEqual(fit.fastKnotPoints, bootstrap.bootstrapMetadata.fastKnotPoints{iBest}, "AbsTol", 1e-12)
            for iDim = 1:3
                testCase.verifyEqual(fit.psiKnotPoints{iDim}, bootstrap.bootstrapMetadata.psiKnotPoints{iBest}{iDim}, "AbsTol", 1e-12)
            end
            testCase.verifyEqual(summary.uCenter, bootstrap.summary.uCenter(:, iBest), "AbsTol", 1e-12)
            testCase.verifyEqual(summary.vCenter, bootstrap.summary.vCenter(:, iBest), "AbsTol", 1e-12)
            testCase.verifyEqual(summary.sigma_n, bootstrap.summary.sigma_n(:, iBest), "AbsTol", 1e-12)
            testCase.verifyEqual(summary.sigma_s, bootstrap.summary.sigma_s(:, iBest), "AbsTol", 1e-12)
            testCase.verifyEqual(summary.zeta, bootstrap.summary.zeta(:, iBest), "AbsTol", 1e-12)
            testCase.verifyEqual(diagnostics.kappa, bootstrap.bootstrapKappa(iBest), "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrap.bestFitKappa, bootstrap.bootstrapKappa(iBest), "AbsTol", 1e-12)
            testCase.verifyTrue(isequaln(bootstrap.bestFitCoherence, bootstrap.bootstrapCoherence(iBest)))
            testCase.verifyEqual(fit.mesoscaleDegreesOfFreedom, bootstrap.mesoscaleDegreesOfFreedom)
            testCase.verifyEqual(bootstrap.bestFit().mesoscaleDegreesOfFreedom, bootstrap.mesoscaleDegreesOfFreedom)
        end

        function mesoscaleDegreesOfFreedomMatchesReconstructedReplicate(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionBootstrapUnitTests.synchronousLinearFieldData();
            bootstrap = GriddedStreamfunctionBootstrap.fromTrajectories(trajectories, nBootstraps=5, randomSeed=13);
            iBootstrap = 3;
            fit = bootstrap.fitForBootstrap(iBootstrap);

            testCase.verifyEqual(bootstrap.mesoscaleDegreesOfFreedom, bootstrap.fullFit.mesoscaleDegreesOfFreedom)
            testCase.verifyEqual(bootstrap.mesoscaleDegreesOfFreedom, fit.mesoscaleDegreesOfFreedom)
        end

        function storedOutputsMatchFullyReconstructedFitsForAllReplicates(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionBootstrapUnitTests.synchronousLinearFieldData();

            bootstrap = GriddedStreamfunctionBootstrap.fromTrajectories(trajectories, nBootstraps=4, randomSeed=11);

            for iBootstrap = 1:bootstrap.nBootstraps
                fit = bootstrap.fitForBootstrap(iBootstrap);
                [summary, diagnostics] = GriddedStreamfunctionBootstrapUnitTests.summaryFromFit(fit, bootstrap.queryTimes);

                testCase.verifyEqual(fit.fastKnotPoints, bootstrap.bootstrapMetadata.fastKnotPoints{iBootstrap}, "AbsTol", 1e-12)
                for iDim = 1:3
                    testCase.verifyEqual(fit.psiKnotPoints{iDim}, ...
                        bootstrap.bootstrapMetadata.psiKnotPoints{iBootstrap}{iDim}, "AbsTol", 1e-12)
                end

                testCase.verifyEqual(summary.uCenter, bootstrap.summary.uCenter(:, iBootstrap), "AbsTol", 1e-12)
                testCase.verifyEqual(summary.vCenter, bootstrap.summary.vCenter(:, iBootstrap), "AbsTol", 1e-12)
                testCase.verifyEqual(summary.sigma_n, bootstrap.summary.sigma_n(:, iBootstrap), "AbsTol", 1e-12)
                testCase.verifyEqual(summary.sigma_s, bootstrap.summary.sigma_s(:, iBootstrap), "AbsTol", 1e-12)
                testCase.verifyEqual(summary.zeta, bootstrap.summary.zeta(:, iBootstrap), "AbsTol", 1e-12)
                testCase.verifyEqual(diagnostics.kappa, bootstrap.bootstrapKappa(iBootstrap), "AbsTol", 1e-12)
                testCase.verifyEqual(fit.mesoscaleDegreesOfFreedom, bootstrap.mesoscaleDegreesOfFreedom)
            end
        end

        function explicitKnotsAreStoredAndReusedDuringReconstruction(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionBootstrapUnitTests.synchronousLinearFieldData();
            referenceFit = GriddedStreamfunction.fromTrajectories(trajectories);

            bootstrap = GriddedStreamfunctionBootstrap.fromTrajectories( ...
                trajectories, ...
                nBootstraps=4, ...
                randomSeed=4, ...
                fastKnotPoints=referenceFit.fastKnotPoints, ...
                psiKnotPoints=referenceFit.psiKnotPoints, ...
                psiS=referenceFit.psiS, ...
                fastS=referenceFit.fastS);

            for iBootstrap = 1:bootstrap.nBootstraps
                testCase.verifyEqual(bootstrap.bootstrapMetadata.fastKnotPoints{iBootstrap}, referenceFit.fastKnotPoints, "AbsTol", 1e-12)
                for iDim = 1:3
                    testCase.verifyEqual(bootstrap.bootstrapMetadata.psiKnotPoints{iBootstrap}{iDim}, ...
                        referenceFit.psiKnotPoints{iDim}, "AbsTol", 1e-12)
                end
            end

            fit = bootstrap.fitForBootstrap(2);
            testCase.verifyEqual(fit.fastKnotPoints, referenceFit.fastKnotPoints, "AbsTol", 1e-12)
            for iDim = 1:3
                testCase.verifyEqual(fit.psiKnotPoints{iDim}, referenceFit.psiKnotPoints{iDim}, "AbsTol", 1e-12)
            end
        end

        function fitForBootstrapMatchesDirectFitWithStoredMetadata(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionBootstrapUnitTests.synchronousLinearFieldData();
            bootstrap = GriddedStreamfunctionBootstrap.fromTrajectories(trajectories, nBootstraps=4, randomSeed=13);
            iBootstrap = 2;
            sampledTrajectories = reshape( ...
                bootstrap.observedTrajectories(bootstrap.bootstrapIndices(iBootstrap, :)), [], 1);

            fit = bootstrap.fitForBootstrap(iBootstrap);
            referenceFit = GriddedStreamfunction.fromTrajectories( ...
                sampledTrajectories, ...
                psiKnotPoints=bootstrap.bootstrapMetadata.psiKnotPoints{iBootstrap}, ...
                psiS=bootstrap.fullFit.psiS, ...
                fastKnotPoints=bootstrap.bootstrapMetadata.fastKnotPoints{iBootstrap}, ...
                fastS=bootstrap.fullFit.fastS, ...
                mesoscaleConstraint=bootstrap.fullFit.mesoscaleConstraint);

            [summary, diagnostics] = GriddedStreamfunctionBootstrapUnitTests.summaryFromFit(fit, bootstrap.queryTimes);
            [referenceSummary, referenceDiagnostics] = GriddedStreamfunctionBootstrapUnitTests.summaryFromFit( ...
                referenceFit, bootstrap.queryTimes);

            testCase.verifyEqual(summary.uCenter, referenceSummary.uCenter, "AbsTol", 1e-12)
            testCase.verifyEqual(summary.vCenter, referenceSummary.vCenter, "AbsTol", 1e-12)
            testCase.verifyEqual(summary.sigma_n, referenceSummary.sigma_n, "AbsTol", 1e-12)
            testCase.verifyEqual(summary.sigma_s, referenceSummary.sigma_s, "AbsTol", 1e-12)
            testCase.verifyEqual(summary.zeta, referenceSummary.zeta, "AbsTol", 1e-12)
            testCase.verifyEqual(diagnostics.kappa, referenceDiagnostics.kappa, "AbsTol", 1e-12)
        end

        function storedBootstrapDiagnosticsMatchDirectReferenceFitForReplicate(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionBootstrapUnitTests.synchronousLinearFieldData();
            bootstrap = GriddedStreamfunctionBootstrap.fromTrajectories(trajectories, nBootstraps=4, randomSeed=13);
            iBootstrap = 2;
            sampledTrajectories = reshape( ...
                bootstrap.observedTrajectories(bootstrap.bootstrapIndices(iBootstrap, :)), [], 1);

            referenceFit = GriddedStreamfunction.fromTrajectories( ...
                sampledTrajectories, ...
                psiKnotPoints=bootstrap.bootstrapMetadata.psiKnotPoints{iBootstrap}, ...
                psiS=bootstrap.fullFit.psiS, ...
                fastKnotPoints=bootstrap.bootstrapMetadata.fastKnotPoints{iBootstrap}, ...
                fastS=bootstrap.fullFit.fastS, ...
                mesoscaleConstraint=bootstrap.fullFit.mesoscaleConstraint);
            [~, referenceDiagnostics] = GriddedStreamfunctionBootstrapUnitTests.summaryFromFit( ...
                referenceFit, bootstrap.queryTimes);

            testCase.verifyEqual(bootstrap.bootstrapKappa(iBootstrap), referenceDiagnostics.kappa, "AbsTol", 1e-12)
            testCase.verifyTrue(isequaln(bootstrap.bootstrapCoherence(iBootstrap), referenceDiagnostics.coherence))
        end

        function fullSummaryMatchesLegacyCenterDiagnostics(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionBootstrapUnitTests.synchronousLinearFieldData();
            bootstrap = GriddedStreamfunctionBootstrap.fromTrajectories(trajectories, nBootstraps=4, randomSeed=7);

            [referenceSummary, referenceDiagnostics] = GriddedStreamfunctionBootstrapUnitTests.summaryFromFit( ...
                bootstrap.fullFit, bootstrap.queryTimes);

            testCase.verifyEqual(bootstrap.fullSummary.uCenter, referenceSummary.uCenter, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrap.fullSummary.vCenter, referenceSummary.vCenter, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrap.fullSummary.sigma_n, referenceSummary.sigma_n, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrap.fullSummary.sigma_s, referenceSummary.sigma_s, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrap.fullSummary.zeta, referenceSummary.zeta, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrap.fullFitKappa, referenceDiagnostics.kappa, "AbsTol", 1e-12)
            testCase.verifyTrue(isequaln(bootstrap.fullFitCoherence, referenceDiagnostics.coherence))
            testCase.verifyEqual(bootstrap.fullFitCoherenceSpectrum.frequency, ...
                referenceDiagnostics.coherenceSpectrum.frequency, "AbsTol", 1e-12)
            testCase.verifyTrue(isequaln(bootstrap.fullFitCoherenceSpectrum.coherence, ...
                referenceDiagnostics.coherenceSpectrum.coherence))
        end

        function fullFitScalarCoherenceUsesOnlyLowerFrequencyHalf(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionBootstrapUnitTests.synchronousLinearFieldData();
            bootstrap = GriddedStreamfunctionBootstrap.fromTrajectories(trajectories, nBootstraps=4, randomSeed=7);

            meanCoherence = bootstrap.fullFitCoherenceSpectrum.coherence;
            lowerFrequencyHalf = meanCoherence(1:ceil(numel(meanCoherence) / 2));
            finiteLowerFrequencyHalf = lowerFrequencyHalf(isfinite(lowerFrequencyHalf));
            expectedCoherence = mean(finiteLowerFrequencyHalf);

            testCase.verifyEqual(bootstrap.fullFitCoherence, expectedCoherence, "AbsTol", 1e-12)

            if numel(meanCoherence) >= 4
                upperFrequencyHalf = meanCoherence((ceil(numel(meanCoherence) / 2) + 1):end);
                finiteUpperFrequencyHalf = upperFrequencyHalf(isfinite(upperFrequencyHalf));
                if ~isempty(finiteUpperFrequencyHalf)
                    modifiedSpectrum = meanCoherence;
                    modifiedSpectrum((ceil(numel(modifiedSpectrum) / 2) + 1):end) = modifiedSpectrum((ceil(numel(modifiedSpectrum) / 2) + 1):end) + 1;
                    modifiedLowerFrequencyHalf = modifiedSpectrum(1:ceil(numel(modifiedSpectrum) / 2));
                    modifiedExpectedCoherence = mean(modifiedLowerFrequencyHalf(isfinite(modifiedLowerFrequencyHalf)));
                    testCase.verifyEqual(modifiedExpectedCoherence, expectedCoherence, "AbsTol", 1e-12)
                end
            end
        end

        function defaultKnotsStillReconstructExactly(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionBootstrapUnitTests.synchronousLinearFieldData();

            bootstrap = GriddedStreamfunctionBootstrap.fromTrajectories(trajectories, nBootstraps=5, randomSeed=12);
            fit = bootstrap.fitForBootstrap(2);
            [summary, diagnostics] = GriddedStreamfunctionBootstrapUnitTests.summaryFromFit(fit, bootstrap.queryTimes);

            testCase.verifyEqual(fit.fastKnotPoints, bootstrap.bootstrapMetadata.fastKnotPoints{2}, "AbsTol", 1e-12)
            for iDim = 1:3
                testCase.verifyEqual(fit.psiKnotPoints{iDim}, bootstrap.bootstrapMetadata.psiKnotPoints{2}{iDim}, "AbsTol", 1e-12)
            end
            testCase.verifyEqual(summary.uCenter, bootstrap.summary.uCenter(:, 2), "AbsTol", 1e-12)
            testCase.verifyEqual(summary.vCenter, bootstrap.summary.vCenter(:, 2), "AbsTol", 1e-12)
            testCase.verifyEqual(summary.sigma_n, bootstrap.summary.sigma_n(:, 2), "AbsTol", 1e-12)
            testCase.verifyEqual(summary.sigma_s, bootstrap.summary.sigma_s(:, 2), "AbsTol", 1e-12)
            testCase.verifyEqual(summary.zeta, bootstrap.summary.zeta(:, 2), "AbsTol", 1e-12)
            testCase.verifyEqual(diagnostics.kappa, bootstrap.bootstrapKappa(2), "AbsTol", 1e-12)
        end

        function defaultQueryTimesUseOriginalCommonOverlap(testCase)
            [trajectories, expectedQueryTimes, overlapBounds] = ...
                GriddedStreamfunctionBootstrapUnitTests.partialSupportTrajectoryData();

            bootstrap = GriddedStreamfunctionBootstrap.fromTrajectories(trajectories, nBootstraps=4, randomSeed=3);

            testCase.verifyEqual(bootstrap.queryTimes, expectedQueryTimes, "AbsTol", 1e-12)
            testCase.verifyGreaterThanOrEqual(min(bootstrap.queryTimes), overlapBounds(1))
            testCase.verifyLessThanOrEqual(max(bootstrap.queryTimes), overlapBounds(2))
            testCase.verifyTrue(all(isfinite(bootstrap.summary.uCenter), "all"))
            testCase.verifyTrue(all(isfinite(bootstrap.summary.sigma_n), "all"))
        end

        function irregularSupportCoherenceUsesCommonUniformOverlapGrid(testCase)
            [trajectories, ~, ~] = GriddedStreamfunctionBootstrapUnitTests.partialSupportTrajectoryData();
            bootstrap = GriddedStreamfunctionBootstrap.fromTrajectories(trajectories, nBootstraps=4, randomSeed=3);
            referenceDiagnostics = GriddedStreamfunctionBootstrapUnitTests.diagnosticsFromFit(bootstrap.fullFit);

            testCase.verifyEqual(bootstrap.fullFitCoherenceSpectrum.frequency, ...
                referenceDiagnostics.coherenceSpectrum.frequency, "AbsTol", 1e-12)
            testCase.verifyTrue(isequaln(bootstrap.fullFitCoherenceSpectrum.coherence, ...
                referenceDiagnostics.coherenceSpectrum.coherence))
        end

        function constrainedBootstrapSkipsDegenerateScoreBlocks(testCase)
            [trajectoriesVorticity, trajectoriesStrain] = ...
                GriddedStreamfunctionBootstrapUnitTests.constrainedTrajectoryData();

            zeroVorticity = GriddedStreamfunctionBootstrap.fromTrajectories( ...
                trajectoriesVorticity, nBootstraps=4, randomSeed=8, mesoscaleConstraint="zeroVorticity");
            zeroStrain = GriddedStreamfunctionBootstrap.fromTrajectories( ...
                trajectoriesStrain, nBootstraps=4, randomSeed=8, mesoscaleConstraint="zeroStrain");

            testCase.verifyEqual(zeroVorticity.scores.zeta, zeros(1, zeroVorticity.nBootstraps), "AbsTol", 1e-12)
            testCase.verifyEqual(zeroVorticity.scores.joint, ...
                zeroVorticity.scores.uv + zeroVorticity.scores.strain, "AbsTol", 1e-12)
            testCase.verifyTrue(all(isfinite(zeroVorticity.scores.joint)))

            testCase.verifyEqual(zeroStrain.scores.strain, zeros(1, zeroStrain.nBootstraps), "AbsTol", 1e-12)
            testCase.verifyEqual(zeroStrain.scores.joint, ...
                zeroStrain.scores.uv + zeroStrain.scores.zeta, "AbsTol", 1e-12)
            testCase.verifyTrue(all(isfinite(zeroStrain.scores.joint)))
        end

        function consensusScoresMatchDirectGaussianReference(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionBootstrapUnitTests.synchronousLinearFieldData();
            bootstrap = GriddedStreamfunctionBootstrap.fromTrajectories(trajectories, nBootstraps=4, randomSeed=5);
            scoreIndices = arrayfun(@(ti) find(bootstrap.queryTimes == ti, 1), bootstrap.scoreTimes);

            expected = GriddedStreamfunctionBootstrapUnitTests.directConsensusScores( ...
                bootstrap.summary, scoreIndices, bootstrap.fullFit.mesoscaleConstraint);

            testCase.verifyEqual(bootstrap.scores.uv, expected.uv, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrap.scores.strain, expected.strain, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrap.scores.zeta, expected.zeta, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrap.scores.joint, expected.joint, "AbsTol", 1e-12)
        end

        function summaryQuantilesReturnExpectedShapesAndMonotoneBands(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionBootstrapUnitTests.synchronousLinearFieldData();

            bootstrap = GriddedStreamfunctionBootstrap.fromTrajectories(trajectories, nBootstraps=5, randomSeed=14);
            quantiles = bootstrap.summaryQuantiles([0.1 0.5 0.9]);

            testCase.verifySize(quantiles.uCenter, [numel(bootstrap.queryTimes), 3])
            testCase.verifySize(quantiles.vCenter, [numel(bootstrap.queryTimes), 3])
            testCase.verifySize(quantiles.sigma_n, [numel(bootstrap.queryTimes), 3])
            testCase.verifySize(quantiles.sigma_s, [numel(bootstrap.queryTimes), 3])
            testCase.verifySize(quantiles.zeta, [numel(bootstrap.queryTimes), 3])
            testCase.verifySize(quantiles.kappa, [1, 3])
            testCase.verifySize(quantiles.coherence, [1, 3])

            testCase.verifyTrue(all(quantiles.uCenter(:, 1) <= quantiles.uCenter(:, 2)))
            testCase.verifyTrue(all(quantiles.uCenter(:, 2) <= quantiles.uCenter(:, 3)))
            testCase.verifyTrue(all(quantiles.sigma_n(:, 1) <= quantiles.sigma_n(:, 2)))
            testCase.verifyTrue(all(quantiles.sigma_n(:, 2) <= quantiles.sigma_n(:, 3)))
            testCase.verifyTrue(all(quantiles.kappa(1) <= quantiles.kappa(2)))
            testCase.verifyTrue(all(quantiles.kappa(2) <= quantiles.kappa(3)))
        end

        function lazyBootstrapCoherenceUsesLowerFrequencyHalf(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionBootstrapUnitTests.synchronousLinearFieldData();
            bootstrap = GriddedStreamfunctionBootstrap.fromTrajectories(trajectories, nBootstraps=4, randomSeed=11);
            iBootstrap = 3;

            fit = bootstrap.fitForBootstrap(iBootstrap);
            referenceDiagnostics = GriddedStreamfunctionBootstrapUnitTests.diagnosticsFromFit(fit);
            meanCoherence = referenceDiagnostics.coherenceSpectrum.coherence;
            lowerFrequencyHalf = meanCoherence(1:ceil(numel(meanCoherence) / 2));
            expectedCoherence = mean(lowerFrequencyHalf(isfinite(lowerFrequencyHalf)));

            testCase.verifyEqual(bootstrap.bootstrapCoherence(iBootstrap), expectedCoherence, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrap.bootstrapCoherence(iBootstrap), referenceDiagnostics.coherence, "AbsTol", 1e-12)
        end

        function scalarDiffusivityDiagnosticIsFiniteAndIncreases(testCase)
            trajectoriesLow = GriddedStreamfunctionBootstrapUnitTests.diffusiveTrajectoryData(0.25, 21);
            trajectoriesHigh = GriddedStreamfunctionBootstrapUnitTests.diffusiveTrajectoryData(2.0, 21);

            bootstrapLow = GriddedStreamfunctionBootstrap.fromTrajectories(trajectoriesLow, nBootstraps=4, randomSeed=5);
            bootstrapHigh = GriddedStreamfunctionBootstrap.fromTrajectories(trajectoriesHigh, nBootstraps=4, randomSeed=5);

            testCase.verifyTrue(isfinite(bootstrapLow.fullFitKappa))
            testCase.verifyTrue(isfinite(bootstrapHigh.fullFitKappa))
            testCase.verifyGreaterThanOrEqual(bootstrapLow.fullFitKappa, 0)
            testCase.verifyGreaterThanOrEqual(bootstrapHigh.fullFitKappa, 0)
            testCase.verifyGreaterThan(bootstrapHigh.fullFitKappa, bootstrapLow.fullFitKappa)
            testCase.verifyGreaterThan(mean(bootstrapHigh.bootstrapKappa), mean(bootstrapLow.bootstrapKappa))
        end

        function scalarDiffusivityDiagnosticMatchesLegacySplineIntegration(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionBootstrapUnitTests.synchronousLinearFieldData();
            bootstrap = GriddedStreamfunctionBootstrap.fromTrajectories(trajectories, nBootstraps=4, randomSeed=5);

            legacyKappa = GriddedStreamfunctionBootstrapUnitTests.legacyKappaFromFit(bootstrap.fullFit);

            testCase.verifyEqual(bootstrap.fullFitKappa, legacyKappa, "AbsTol", 1e-12)
        end

        function scalarDiffusivityDiagnosticMatchesDirectTotalDisplacementFormula(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionBootstrapUnitTests.synchronousLinearFieldData();
            bootstrap = GriddedStreamfunctionBootstrap.fromTrajectories(trajectories, nBootstraps=4, randomSeed=5);

            directKappa = GriddedStreamfunctionBootstrapUnitTests.totalDisplacementKappaFromSubmesoscaleTrajectories(bootstrap.fullFit);

            testCase.verifyEqual(bootstrap.fullFitKappa, directKappa, "AbsTol", 1e-12)
        end

        function missingCoherenceBackendWarnsOnceAndLeavesCoherenceUnavailable(testCase)
            if ~GriddedStreamfunctionBootstrapUnitTests.hasCoherenceBackend()
                return
            end

            sleptapPath = fileparts(which('sleptap'));
            mspecPath = fileparts(which('mspec'));
            removedPaths = unique({sleptapPath, mspecPath}, 'stable');
            clear sleptap mspec
            for iPath = 1:numel(removedPaths)
                rmpath(removedPaths{iPath});
            end
            cleanup = onCleanup(@() GriddedStreamfunctionBootstrapUnitTests.restorePaths(removedPaths));

            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionBootstrapUnitTests.synchronousLinearFieldData();
            lastwarn('');
            bootstrap = GriddedStreamfunctionBootstrap.fromTrajectories(trajectories, nBootstraps=4, randomSeed=5);
            [~, warningId] = lastwarn;
            testCase.verifyEqual(string(warningId), "GriddedStreamfunctionBootstrap:MissingCoherenceBackend")

            lastwarn('');
            bootstrap.bootstrapCoherence;
            [~, warningId] = lastwarn;
            testCase.verifyEqual(string(warningId), "")
            testCase.verifyTrue(all(isnan(bootstrap.bootstrapCoherence)))
            testCase.verifyTrue(isnan(bootstrap.fullFitCoherence))
            testCase.verifyTrue(isnan(bootstrap.bestFitCoherence))
            testCase.verifyEmpty(bootstrap.fullFitCoherenceSpectrum.frequency)
            testCase.verifyEmpty(bootstrap.bestFitCoherenceSpectrum.frequency)
        end
    end

    methods (Static)
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
            trajectories = GriddedStreamfunctionBootstrapUnitTests.trajectorySplinesFromMatrices(t, x, y);
        end

        function [trajectories, expectedQueryTimes, overlapBounds] = partialSupportTrajectoryData()
            sigma = 3.0e-6;
            theta = pi/10;
            zeta = -1.0e-6;
            model = LinearVelocityField(sigma=sigma, theta=theta, zeta=zeta, u0=0.04, v0=-0.02);
            integrator = AdvectionDiffusionIntegrator(model, 0);
            [x0, y0] = ndgrid([-400 0 400] + 800, [-300 0 300] - 500);
            [tFine, xFine, yFine] = integrator.particleTrajectories(x0(:), y0(:), 10*3600, 300);

            startIndices = [1; 4; 2; 6; 3];
            endIndices = [numel(tFine) - 2; numel(tFine); numel(tFine) - 4; numel(tFine) - 1; numel(tFine) - 3];
            nDrifters = numel(startIndices);
            trajectories = TrajectorySpline.empty(0, 1);
            allTimes = cell(nDrifters, 1);

            for iDrifter = 1:nDrifters
                indices = startIndices(iDrifter):3:endIndices(iDrifter);
                allTimes{iDrifter} = tFine(indices);
                trajectories(end + 1, 1) = TrajectorySpline.fromData(tFine(indices), xFine(indices, iDrifter), yFine(indices, iDrifter), S=3); %#ok<AGROW>
            end

            overlapBounds = [max(cellfun(@(ti) ti(1), allTimes)), min(cellfun(@(ti) ti(end), allTimes))];
            pooledTimes = unique(vertcat(allTimes{:}), "sorted");
            expectedQueryTimes = pooledTimes(pooledTimes >= overlapBounds(1) & pooledTimes <= overlapBounds(2));
        end

        function [trajectoriesZeroVorticity, trajectoriesZeroStrain] = constrainedTrajectoryData()
            integratorVorticity = AdvectionDiffusionIntegrator(LinearVelocityField(sigma=4.0e-6, theta=pi/8, zeta=0, u0=0.06, v0=-0.03), 0);
            integratorStrain = AdvectionDiffusionIntegrator(LinearVelocityField(sigma=0, theta=pi/8, zeta=-2.0e-6, u0=0.06, v0=-0.03), 0);
            [x0, y0] = ndgrid([-600 0 600] + 1200, [-400 0 400] - 800);
            [t, xVorticity, yVorticity] = integratorVorticity.particleTrajectories(x0(:), y0(:), 12*3600, 900);
            [~, xStrain, yStrain] = integratorStrain.particleTrajectories(x0(:), y0(:), 12*3600, 900);

            trajectoriesZeroVorticity = GriddedStreamfunctionBootstrapUnitTests.trajectorySplinesFromMatrices(t, xVorticity, yVorticity);
            trajectoriesZeroStrain = GriddedStreamfunctionBootstrapUnitTests.trajectorySplinesFromMatrices(t, xStrain, yStrain);
        end

        function trajectories = diffusiveTrajectoryData(kappa, rngSeed)
            sigma = 3.0e-6;
            theta = pi/12;
            zeta = -1.0e-6;
            model = LinearVelocityField(sigma=sigma, theta=theta, zeta=zeta, u0=0.03, v0=-0.02);
            integrator = AdvectionDiffusionIntegrator(model, kappa);
            integrator.stepSize = 900;

            [x0, y0] = ndgrid([-500 0 500] + 1000, [-350 0 350] - 600);
            rng(rngSeed);
            [t, x, y] = integrator.particleTrajectories(x0(:), y0(:), 18*3600, 1800);
            trajectories = GriddedStreamfunctionBootstrapUnitTests.trajectorySplinesFromMatrices(t, x, y);
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

        function [summary, diagnostics] = summaryFromFit(fit, queryTimes)
            mx = reshape(fit.centerOfMassTrajectory.x(queryTimes), [], 1);
            my = reshape(fit.centerOfMassTrajectory.y(queryTimes), [], 1);

            summary = struct( ...
                "uCenter", reshape(fit.uMesoscale(queryTimes, mx, my), [], 1), ...
                "vCenter", reshape(fit.vMesoscale(queryTimes, mx, my), [], 1), ...
                "sigma_n", reshape(fit.sigma_n(queryTimes, mx, my), [], 1), ...
                "sigma_s", reshape(fit.sigma_s(queryTimes, mx, my), [], 1), ...
                "zeta", reshape(fit.zeta(queryTimes, mx, my), [], 1));
            diagnostics = GriddedStreamfunctionBootstrapUnitTests.diagnosticsFromFit(fit);
        end

        function kappaEstimate = kappaFromFit(fit)
            nTrajectories = numel(fit.decomposition.fixedFrame.submesoscale);
            kappaValues = zeros(nTrajectories, 1);

            for iTrajectory = 1:nTrajectories
                trajectory = fit.decomposition.fixedFrame.submesoscale(iTrajectory);
                ti = trajectory.t;
                duration = ti(end) - ti(1);
                kappaValues(iTrajectory) = (trajectory.x(ti(end)).^2 + trajectory.y(ti(end)).^2) / (4 * duration);
            end

            kappaEstimate = mean(kappaValues);
        end

        function diagnostics = diagnosticsFromFit(fit)
            diagnostics = struct( ...
                "kappa", GriddedStreamfunctionBootstrapUnitTests.kappaFromFit(fit), ...
                "coherence", NaN, ...
                "coherenceSpectrum", struct("frequency", zeros(0, 1), "coherence", zeros(0, 1)));

            if ~GriddedStreamfunctionBootstrapUnitTests.hasCoherenceBackend()
                return
            end

            trajectories = reshape(fit.observedTrajectories, [], 1);
            tCell = arrayfun(@(trajectory) reshape(trajectory.t, [], 1), trajectories, 'UniformOutput', false);
            tCoherence = GriddedStreamfunctionBootstrapUnitTests.coherenceEvaluationTimes(tCell);
            if numel(tCoherence) < 2
                return
            end

            centeredFrame = fit.decomposition.centeredFrame;
            nTrajectories = numel(centeredFrame.mesoscale);
            cvMesoscale = complex(zeros(numel(tCoherence), nTrajectories));
            cvSubmesoscale = complex(zeros(numel(tCoherence), nTrajectories));

            for iTrajectory = 1:nTrajectories
                mesoscale = centeredFrame.mesoscale(iTrajectory);
                submesoscale = centeredFrame.submesoscale(iTrajectory);
                cvMesoscale(:, iTrajectory) = reshape(mesoscale.u(tCoherence), [], 1) + 1i * reshape(mesoscale.v(tCoherence), [], 1);
                cvSubmesoscale(:, iTrajectory) = reshape(submesoscale.u(tCoherence), [], 1) + 1i * reshape(submesoscale.v(tCoherence), [], 1);
            end

            [psi, ~] = sleptap(size(cvMesoscale, 1));
            dtCoherence = tCoherence(2) - tCoherence(1);
            [frequency, sxx, syy, sxy] = mspec(dtCoherence, cvMesoscale, cvSubmesoscale, psi, 'cyclic');
            gamma = abs(sxy).^2 ./ (sxx .* syy);
            meanCoherence = GriddedStreamfunctionBootstrapUnitTests.meanOverFinite(gamma);
            frequency = reshape(frequency, [], 1);
            meanCoherence = reshape(meanCoherence, [], 1);
            lowerFrequencyHalf = meanCoherence(1:ceil(numel(meanCoherence) / 2));
            finiteCoherence = isfinite(lowerFrequencyHalf);
            if ~any(finiteCoherence)
                return
            end

            diagnostics.coherence = mean(lowerFrequencyHalf(finiteCoherence));
            diagnostics.coherenceSpectrum = struct( ...
                "frequency", frequency, ...
                "coherence", meanCoherence);
        end

        function tf = hasCoherenceBackend()
            tf = exist('sleptap', 'file') ~= 0 && exist('mspec', 'file') ~= 0;
        end

        function tCoherence = coherenceEvaluationTimes(tCell)
            nTrajectories = numel(tCell);
            tStart = -Inf;
            tEnd = Inf;
            dtValues = zeros(nTrajectories, 1);

            for iTrajectory = 1:nTrajectories
                ti = reshape(tCell{iTrajectory}, [], 1);
                if numel(ti) < 2
                    tCoherence = zeros(0, 1);
                    return
                end

                dtValues(iTrajectory) = median(diff(ti));
                if ~isfinite(dtValues(iTrajectory)) || dtValues(iTrajectory) <= 0
                    tCoherence = zeros(0, 1);
                    return
                end

                tStart = max(tStart, ti(1));
                tEnd = min(tEnd, ti(end));
            end

            if ~(isfinite(tStart) && isfinite(tEnd)) || tEnd <= tStart
                tCoherence = zeros(0, 1);
                return
            end

            dtCoherence = max(dtValues);
            nStep = floor((tEnd - tStart) / dtCoherence);
            if nStep < 1
                tCoherence = zeros(0, 1);
                return
            end

            tCoherence = tStart + (0:nStep).' * dtCoherence;
        end

        function meanValues = meanOverFinite(values)
            meanValues = zeros(size(values, 1), 1);
            for iRow = 1:size(values, 1)
                finiteValues = values(iRow, isfinite(values(iRow, :)));
                if isempty(finiteValues)
                    meanValues(iRow) = NaN;
                else
                    meanValues(iRow) = mean(finiteValues);
                end
            end
        end

        function kappaEstimate = legacyKappaFromFit(fit)
            trajectories = reshape(fit.observedTrajectories, [], 1);
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

            allT = vertcat(tCell{:});
            allX = vertcat(xCell{:});
            allY = vertcat(yCell{:});
            allXDot = vertcat(xDotCell{:});
            allYDot = vertcat(yDotCell{:});
            mxDotAll = reshape(fit.centerOfMassTrajectory.u(allT), [], 1);
            myDotAll = reshape(fit.centerOfMassTrajectory.v(allT), [], 1);
            uMesoscaleObserved = reshape(fit.uMesoscale(allT, allX, allY), [], 1);
            vMesoscaleObserved = reshape(fit.vMesoscale(allT, allX, allY), [], 1);
            B = reshape(BSpline.matrixForDataPoints(allT, knotPoints=fit.fastKnotPoints, S=fit.fastS), numel(allT), []);
            uMesoscaleComObserved = B * (B \ uMesoscaleObserved);
            vMesoscaleComObserved = B * (B \ vMesoscaleObserved);
            uBackgroundObserved = mxDotAll - uMesoscaleComObserved;
            vBackgroundObserved = myDotAll - vMesoscaleComObserved;
            uSubmesoscaleObserved = allXDot - uBackgroundObserved - uMesoscaleObserved;
            vSubmesoscaleObserved = allYDot - vBackgroundObserved - vMesoscaleObserved;

            kappaValues = zeros(nTrajectories, 1);
            sampleStartIndex = 1;

            for iTrajectory = 1:nTrajectories
                ti = tCell{iTrajectory};
                duration = ti(end) - ti(1);
                nSamples = numel(ti);
                sampleIndices = sampleStartIndex:(sampleStartIndex + nSamples - 1);
                componentS = min(3, nSamples - 1);
                uSpline = InterpolatingSpline.fromGriddedValues(ti, uSubmesoscaleObserved(sampleIndices), S=componentS);
                vSpline = InterpolatingSpline.fromGriddedValues(ti, vSubmesoscaleObserved(sampleIndices), S=componentS);
                xEnd = cumsum(uSpline);
                yEnd = cumsum(vSpline);
                kappaValues(iTrajectory) = (xEnd(ti(end)).^2 + yEnd(ti(end)).^2) / (4 * duration);
                sampleStartIndex = sampleStartIndex + nSamples;
            end

            kappaEstimate = mean(kappaValues);
        end

        function kappaEstimate = totalDisplacementKappaFromSubmesoscaleTrajectories(fit)
            trajectories = reshape(fit.decomposition.fixedFrame.submesoscale, [], 1);
            nTrajectories = numel(trajectories);
            kappaValues = zeros(nTrajectories, 1);

            for iTrajectory = 1:nTrajectories
                trajectory = trajectories(iTrajectory);
                ti = trajectory.t;
                duration = ti(end) - ti(1);
                componentS = min(3, numel(ti) - 1);
                tKnot = BSpline.knotPointsForDataPoints(ti, S=componentS);
                totalDisplacementWeights = BSpline.integralMatrixForDataPoints(ti, ti(end), knotPoints=tKnot, S=componentS);
                xTotal = totalDisplacementWeights * reshape(trajectory.u(ti), [], 1);
                yTotal = totalDisplacementWeights * reshape(trajectory.v(ti), [], 1);
                kappaValues(iTrajectory) = (xTotal.^2 + yTotal.^2) / (4 * duration);
            end

            kappaEstimate = mean(kappaValues);
        end

        function scores = directConsensusScores(summary, scoreIndices, mesoscaleConstraint)
            nBootstraps = size(summary.uCenter, 2);
            scores = struct( ...
                "uv", zeros(1, nBootstraps), ...
                "strain", zeros(1, nBootstraps), ...
                "zeta", zeros(1, nBootstraps), ...
                "joint", zeros(1, nBootstraps));

            for iScore = reshape(scoreIndices, 1, [])
                scores.uv = scores.uv + GriddedStreamfunctionBootstrapUnitTests.directScore2DBlock( ...
                    summary.uCenter(iScore, :), summary.vCenter(iScore, :));

                if mesoscaleConstraint ~= "zeroStrain"
                    scores.strain = scores.strain + GriddedStreamfunctionBootstrapUnitTests.directScore2DBlock( ...
                        summary.sigma_n(iScore, :), summary.sigma_s(iScore, :));
                end

                if mesoscaleConstraint ~= "zeroVorticity"
                    scores.zeta = scores.zeta + GriddedStreamfunctionBootstrapUnitTests.directScore1DBlock( ...
                        summary.zeta(iScore, :));
                end
            end

            scores.joint = scores.uv + scores.strain + scores.zeta;
        end

        function score = directScore2DBlock(xValues, yValues)
            xValues = reshape(xValues, [], 1);
            yValues = reshape(yValues, [], 1);
            score = zeros(1, numel(xValues));
            [isConstantX, isConstantY] = GriddedStreamfunctionBootstrapUnitTests.constantFlags(xValues, yValues);

            if isConstantX && isConstantY
                return
            end

            if isConstantX
                score = GriddedStreamfunctionBootstrapUnitTests.directScore1DBlock(yValues);
                return
            end

            if isConstantY
                score = GriddedStreamfunctionBootstrapUnitTests.directScore1DBlock(xValues);
                return
            end

            model = KernelDensityEstimate.fromData([xValues, yValues]);
            pointDensity = GriddedStreamfunctionBootstrapUnitTests.directGaussianDensity2D(model.data, model.bandwidth, [xValues, yValues]);
            score = log10(GriddedStreamfunctionBootstrapUnitTests.safeDensity(pointDensity)).';
        end

        function score = directScore1DBlock(values)
            values = reshape(values, [], 1);
            score = zeros(1, numel(values));

            if GriddedStreamfunctionBootstrapUnitTests.isConstant(values)
                return
            end

            model = KernelDensityEstimate.fromData(values);
            pointDensity = GriddedStreamfunctionBootstrapUnitTests.directGaussianDensity1D(model.data, model.bandwidth, values);
            score = log10(GriddedStreamfunctionBootstrapUnitTests.safeDensity(pointDensity)).';
        end

        function density = safeDensity(values)
            density = reshape(values, [], 1);
            density(~isfinite(density) | density <= 0) = realmin("double");
        end

        function [isConstantX, isConstantY] = constantFlags(xValues, yValues)
            isConstantX = GriddedStreamfunctionBootstrapUnitTests.isConstant(xValues);
            isConstantY = GriddedStreamfunctionBootstrapUnitTests.isConstant(yValues);
        end

        function tf = isConstant(values)
            values = reshape(values, [], 1);
            scale = max([1; abs(values)]);
            tf = (max(values) - min(values)) <= 1e-12 * scale;
        end

        function density = directGaussianDensity1D(data, bandwidth, queryPoints)
            queryPoints = reshape(queryPoints, [], 1);
            density = mean(exp(-0.5 * ((queryPoints - data.')/bandwidth).^2), 2)/(bandwidth * sqrt(2*pi));
        end

        function density = directGaussianDensity2D(data, bandwidth, queryPoints)
            density = mean(exp(-0.5 * (((queryPoints(:, 1) - data(:, 1).')/bandwidth(1)).^2 + ((queryPoints(:, 2) - data(:, 2).')/bandwidth(2)).^2)), 2) ...
                /(2*pi * bandwidth(1) * bandwidth(2));
        end

        function deleteFileIfExists(path)
            if isfile(path)
                delete(path);
            end
        end

        function deleteDirectoryIfExists(path)
            if exist(path, 'dir') ~= 0
                rmdir(path, 's');
            end
        end

        function restorePaths(paths)
            for iPath = 1:numel(paths)
                addpath(paths{iPath});
            end
            clear sleptap mspec
        end
    end
end
