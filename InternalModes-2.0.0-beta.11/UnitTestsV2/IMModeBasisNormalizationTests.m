classdef IMModeBasisNormalizationTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function configurePath(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,"UnitTestsV2","TestSupport")));
        end
    end
    methods (Test)
        function waveRulesAgreeWithFreshConstructionAcrossBackends(testCase)
            bg = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid"); z = linspace(-1000,0,31).';
            for solver = {IMSolverSpectral(nEVP=48),IMSolverFiniteDifference(z=linspace(-1000,0,81).'),IMSolverAnalytical()}
                original = createIMTestBasis("wave",bg,solver=solver{1});
                oldF = original.F(z);
                for name = ["rmsF","rmsG","maxAbsF","maxAbsG","waveEnergy"]
                    changed = original.withNormalization(name); fresh = createIMTestBasis("wave",bg,solver=solver{1},normalization=name);
                    testCase.verifyEqual(changed.F(z),fresh.F(z),AbsTol=3e-12,RelTol=3e-12);
                    testCase.verifyEqual(changed.G(z),fresh.G(z),AbsTol=3e-12,RelTol=3e-12);
                    testCase.verifyEqual(changed.F(z),bg.D*changed.Gz(z),AbsTol=3e-13);
                    testCase.verifyEqual(changed.Fz(z),bg.D*changed.Gzz(z),AbsTol=3e-14);
                    verifyPreserved(testCase,original,changed);
                end
                testCase.verifyEqual(original.F(z),oldF,AbsTol=0);
            end
        end
        function signedEnergyRulesPreserveMetricSignatures(testCase)
            bg = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface"); energy = IMGeneralizedEnergyBoundary(g0=-.04,gd=.03);
            b = createIMTestBasis("MDA",bg,nModes=4,energyBoundary=energy,normalization="meanDensityAnomalyEnergy");
            testCase.verifyTrue(any(b.normSignatures < 0));
            normalized = b.withNormalization("rmsG"); testCase.verifyEqual(normalized.normSignatures,ones(1,b.nColumns));
            restored = normalized.withNormalization("meanDensityAnomalyEnergy");
            testCase.verifyEqual(restored.normSignatures,b.normSignatures);
            custom = b.addNormalization("customScale",@(~,~) 2).withNormalization("customScale");
            testCase.verifyEqual(custom.normSignatures,ones(1,b.nColumns));
            testCase.verifyTrue(any(diag(custom.gramMatrix(variable="G")) < 0));
            testCase.verifyEqual(restored.G([-950;-510;-40]),b.G([-950;-510;-40]),AbsTol=2e-13);
        end
        function remainingEVPFamiliesAgreeWithTheirFactoryRules(testCase)
            bg = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface"); energy = IMGeneralizedEnergyBoundary(g0=.03,gd=.01);
            for family = ["inertial","APV","hydrostaticG"]
                b = createIMTestBasis(family,bg,energyBoundary=energy);
                for name = ["rmsF","maxAbsG"]
                    c = b.withNormalization(name); fresh = createIMTestBasis(family,bg,energyBoundary=energy,normalization=name);
                    testCase.verifyEqual(c.F([-900;-500;-100]),fresh.F([-900;-500;-100]),AbsTol=3e-12,RelTol=3e-12);
                    verifyPreserved(testCase,b,c);
                end
            end
        end
        function generalizedEnstrophyRulesReuseTheSameSelectedModes(testCase)
            bg = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            b = createIMTestBasis("enstrophy",bg,nModes=4); selected = b.selectColumns([4 2]);
            changed = selected.withNormalization("rmsF"); expected = b.withNormalization("rmsF").selectColumns([4 2]);
            testCase.verifyEqual(changed.F([-900;-500;-100]),expected.F([-900;-500;-100]),AbsTol=2e-12,RelTol=2e-12);
            verifyPreserved(testCase,selected,changed);
        end
        function selectedHydrostaticModesDoNotNormalizeDiscardedNullColumns(testCase)
            bg = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            for solver = {IMSolverSpectral(nEVP=48),IMSolverAnalytical()}
                full = createIMTestBasis("hydrostaticF",bg,nModes=4,solver=solver{1});
                invalid = full.withNormalization("rmsG");
                if isa(full,"IMInternalModesBasis"), identifier = "IMInternalModesBasis:NullNormalization"; else, identifier = "IMAnalyticalInternalModesBasis:NullNormalization"; end
                testCase.verifyError(@() invalid.normalizationFactors(),identifier);
                selected = full.selectColumns([4 2]).withNormalization("rmsG");
                [z,w] = legpts(100,bg.zDomain); w = w(:);
                testCase.verifyEqual(sum(w.*selected.G(z).^2,1)/bg.D,ones(1,2),AbsTol=3e-12);
            end
        end
        function neutralLabelsScaleEveryChannelAndRetainDDIdentities(testCase)
            bg = IMBackground(N2=@(z) zeros(size(z)),zDomain=[-1000 0],surfaceCondition="freeSurface");
            full = IMNeutralLabelModesBasis(bg,neutralIntervals=[-1000 -800;-200 0],nModes=[3 2]);
            b = full.selectColumns([5 1 4]); c = b.withNormalization("maxAbsG");
            z = [-950;-850;-150;-50]; ratio = b.normalizationFactors()./c.normalizationFactors();
            testCase.verifyEqual(c.G(z),b.G(z).*ratio,AbsTol=2e-14);
            testCase.verifyEqual(c.APV(z),b.APV(z).*ratio,AbsTol=2e-14);
            testCase.verifyEqual(c.G(b.neutralIntervals(:)),zeros(4,3),AbsTol=0);
            replaced = full.addNormalization("customScale",@(~,j) j+1).selectColumns([5 1 4]).withNormalization("customScale");
            testCase.verifyEqual(replaced.G(z),b.G(z)./[6 2 5],AbsTol=2e-14);
            testCase.verifyEqual(full.normalization,"rmsG");
        end
        function highSelectedHarmonicsKeepTheirNormalizationResolution(testCase)
            bg = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            b = createIMTestBasis("wave",bg,nModes=160,solver=IMSolverAnalytical(),normalization="rmsF");
            selected = b.selectColumns([160 137]).withNormalization("rmsG");
            full = b.withNormalization("rmsG"); expected = full.selectColumns([160 137]);
            testCase.verifyEqual(selected.normalizationFactors(),expected.normalizationFactors(),RelTol=5e-13);
            [z,w] = legpts(700,bg.zDomain); w = w(:);
            testCase.verifyEqual(sum(w.*selected.G(z).^2,1)/bg.D,ones(1,2),AbsTol=5e-13);
            critical = IMBackground(N2=@(z) bg.f0^2*ones(size(z)),zDomain=bg.zDomain,surfaceCondition="freeSurface");
            kernel = IMZeroDepthWaveModesBasis(critical,1e-3,criticalIntervals=[-900 -700],nModesPerInterval=160);
            a = kernel.selectColumns([137 100]).withNormalization("rmsG"); c = kernel.withNormalization("rmsG").selectColumns([137 100]);
            testCase.verifyEqual(a.G([-899;-851;-799;-741;-701]),c.G([-899;-851;-799;-741;-701]),AbsTol=2e-12,RelTol=5e-13);
            replaced = kernel.addNormalization("customScale",@(~,j) j+1).selectColumns([137 100]).withNormalization("customScale");
            testCase.verifyEqual(replaced.G([-850;-750]),a.G([-850;-750])./[138 101],AbsTol=2e-14);
        end
        function fixedFrequencyRootsAndEmptyInventoriesKeepTheirContracts(testCase)
            bg = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            b = createIMTestBasis("frequency",bg,omega=.5*abs(bg.f0),solver=IMSolverAnalytical()); empty = b.selectColumns([]);
            for name = reshape(empty.normalizationNames(),1,[])
                changed = empty.withNormalization(name);
                testCase.verifyEmpty(changed.normalizationFactors());
                testCase.verifyEmpty(changed.normSignatures);
                testCase.verifySize(changed.F([-1000;0]),[2 0]);
            end
            testCase.verifyEqual(empty.kSquared,zeros(1,0));
        end
        function endpointResponsesReuseTheExistingBoundaryValueSolve(testCase)
            bg = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            b = IMSolverAnalytical().solveGeostrophicZeroAPVModes(IMGeostrophicZeroAPVModes.atWavenumber(bg,1e-3,endpoints=["surface","bottom"]));
            changed = b.withNormalization("maxAbsG"); z = linspace(-1000,0,31).';
            ratio = reshape(b.normalizationFactors()./changed.normalizationFactors(),1,[]);
            testCase.verifyEqual(changed.G(z),b.G(z).*ratio,AbsTol=3e-12,RelTol=3e-12);
            testCase.verifyEqual(changed.constructionRecord.requested,b.constructionRecord.requested);
            testCase.verifyEqual(changed.boundaryResponseMatrix,b.boundaryResponseMatrix.*ratio,AbsTol=3e-12,RelTol=3e-12);
        end
        function sameNameAndInvalidRequestsAvoidProfileEvaluation(testCase)
            calls = 0;
            p = IMInternalModes(N2=@(z) ones(size(z)),zDomain=[-1 0],formulation="G",p=1,q=0,r=1,FfromGz=@diagnostic);
            b = IMSolverSpectral(nEVP=32).solveEVP(p,nModes=2); calls = 0;
            changed = b.withNormalization(b.normalization); changed.capabilities; changed.scientificSpecification; changed.constructionRecord;
            testCase.verifyEqual(calls,0);
            testCase.verifyError(@() b.withNormalization("unknownRule"),"IMBasisSet:UnsupportedNormalization");
            testCase.verifyEqual(calls,0);
            function values = diagnostic(~,dG,~,~)
                calls = calls+1; values = dG.^2+1;
            end
        end
        function nonfiniteMagnitudeIsRejectedWithoutChangingTheInput(testCase)
            bg = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            b = createIMTestBasis("wave",bg).addNormalization("bad",@(~,~) Inf); c = b.withNormalization("bad");
            testCase.verifyError(@() c.normalizationFactors(),"IMBasisSet:NullNormalization");
            testCase.verifyEqual(b.normalization,"waveEnergy");
            testCase.verifyTrue(all(isfinite(b.F(-500))));
            bad = b.addNormalization("negative",@(~,~) -2).withNormalization("negative");
            testCase.verifyError(@() bad.normalizationFactors(),"IMBasisSet:NullNormalization");
        end
        function normalizationDoesNotEnterEVPOrWeightFittingConstruction(testCase)
            bg = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface"); b = createIMTestBasis("wave",bg);
            cleanup = onCleanup(@() profile("off")); profile("clear"); profile("on");
            c = b.withNormalization("rmsG"); c.F([-900;-100]); profile("off");
            info = profile("info"); names = string({info.FunctionTable.FunctionName});
            forbidden = ["solveAssembledEVP","solveEVP","solveGeostrophicZeroAPVModes","quadratureWeightsForPoints","IMQuadratureWeightFit"];
            testCase.verifyFalse(any(contains(names,forbidden)));
            verifyPreserved(testCase,b,c); clear cleanup
        end
    end
end
function verifyPreserved(testCase,original,changed)
testCase.verifyEqual(changed.columnLabels,original.columnLabels);
testCase.verifyEqual(changed.modeInfo,original.modeInfo);
testCase.verifyEqual(changed.h,original.h);
testCase.verifyEqual(changed.k,original.k);
testCase.verifyEqual(changed.kSquared,original.kSquared);
testCase.verifyEqual(changed.constructionRecord,original.constructionRecord);
testCase.verifyEqual(changed.scientificSpecification.requested,original.scientificSpecification.requested);
end
