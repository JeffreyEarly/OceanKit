classdef IMNormalizationCatalogueTests < matlab.unittest.TestCase
    % Verify dimensionless defaults and the closed v2 normalization catalogue.

    properties (TestParameter)
        customRuleFamily = {"scalar","numerical","analytical","neutral","zeroDepth"}
    end

    methods (TestClassSetup)
        function configurePath(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,"UnitTestsV2","TestSupport")));
        end
    end

    methods (Test)
        function customRuleDivisorsAreValidatedBeforeConversion(testCase,customRuleFamily)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            switch customRuleFamily
                case "scalar"
                    basis = IMSolverSpectral(nEVP=24).solveEVP(IMEigenvalueProblem(zDomain=[-1 0],p=1,q=0,r=1),nModes=2);
                case "numerical"
                    basis = createIMTestBasis("wave",background,nModes=2,solver=IMSolverSpectral(nEVP=24));
                case "analytical"
                    basis = createIMTestBasis("wave",background,nModes=2,solver=IMSolverAnalytical());
                case "neutral"
                    basis = IMNeutralLabelModesBasis([-1000 0],zDomain=background.zDomain,nModesPerInterval=2);
                case "zeroDepth"
                    critical = IMBackground(N2=@(z) background.f0^2*ones(size(z)),zDomain=background.zDomain,surfaceCondition="rigidLid");
                    basis = IMZeroDepthWaveModesBasis(critical,1e-3,criticalIntervals=[-1000 0],nModesPerInterval=2);
            end
            prefix = string(class(basis));
            if isa(basis,"IMBasisSet"), prefix = "IMBasisSet"; end
            for item = {'x',"2",true,struct(),[],[1 2]}
                value = item{1};
                custom = basis.addNormalization("customScale",@(~,~) value).withNormalization("customScale");
                testCase.verifyError(@() custom.normalizationFactors(),prefix+":InvalidNormalizationRule");
            end
            for value = [2+1i,0,-1,NaN,Inf,-Inf]
                custom = basis.addNormalization("customScale",@(~,~) value).withNormalization("customScale");
                testCase.verifyError(@() custom.normalizationFactors(),prefix+":NullNormalization");
            end
            for item = {2,single(2),uint8(2)}
                value = item{1};
                custom = basis.addNormalization("customScale",@(~,~) value).withNormalization("customScale");
                testCase.verifyEqual(custom.normalizationFactors(),[2 2]);
            end
        end

        function scalarRmsUsesDepthAndPreservesTheScientificMetric(testCase)
            evp = IMEigenvalueProblem(zDomain=[-3 0],p=1,q=0,r=2);
            for solver = {IMSolverSpectral(nEVP=40),IMSolverFiniteDifference(z=linspace(-3,0,121).')}
                basis = solver{1}.solveEVP(evp,nModes=3);
                if isa(basis.solver,"IMSolverFiniteDifference")
                    knots = sort(basis.solver.zNative);
                    [nodes,rule] = legpts(4,[0 1]);
                    z = reshape(knots(1:end-1).'+nodes(:).*diff(knots).',[],1);
                    weights = reshape(rule(:).*diff(knots).',[],1);
                else
                    [z,weights] = legpts(96,basis.zDomain);
                    weights = weights(:);
                end
                testCase.verifyEqual(basis.normalization,"rms");
                testCase.verifyEqual(sum(weights.*basis.u(z).^2,1)/3,ones(1,3),AbsTol=3e-13);
                testCase.verifyEqual(diag(basis.gramMatrix()),6*ones(3,1),AbsTol=2e-12);
                testCase.verifyEqual(basis.normSignatures,ones(1,3));
            end
        end

        function scalarRmsResolvesTheSquareOfTheHighestChebyshevProfile(testCase)
            degree = 23;
            evp = IMEigenvalueProblem(zDomain=[-3 0],p=1,q=0,r=1);
            solver = IMSolverSpectral(nEVP=degree+1).configuredForEVP(evp);
            coefficients = zeros(degree+1,1);
            coefficients(end) = 1;
            % The exactly represented T_degree profile isolates quadrature
            % accuracy from the convergence of a solved eigenfunction.
            basis = IMBasisSet(solver=solver,evp=evp,nativeModes=coefficients,eigenvalues=1,modeNumber=1);
            exactRms = sqrt((1-1/(4*degree^2-1))/2);
            testCase.verifyEqual(basis.normalizationFactors(),exactRms,AbsTol=3e-14);
            measuredRms = sqrt(integral(@(z) abs(basis.u(z)).^2,-3,0,ArrayValued=true,AbsTol=1e-12,RelTol=1e-12)/3);
            testCase.verifyEqual(measuredRms,1,AbsTol=3e-12);
        end

        function scalarRmsNearCapacityAgreesWithIndependentIntegration(testCase)
            evp = IMEigenvalueProblem(zDomain=[-1 0],p=1,q=0,r=1,surfaceBoundary=IMBoundaryCondition.dirichlet(),bottomBoundary=IMBoundaryCondition.dirichlet());
            basis = IMSolverSpectral(nEVP=24).solveEVP(evp,nModes=22);
            % Adaptive physical-coordinate integration is independent of
            % the fixed Gauss rule used to prepare normalization divisors.
            measuredRms = sqrt(integral(@(z) abs(basis.u(z)).^2,-1,0,ArrayValued=true,AbsTol=1e-12,RelTol=1e-12));
            testCase.verifyEqual(measuredRms,ones(1,22),AbsTol=3e-12);
        end

        function pairedRmsRetainsPhysicalDepthMeasureInStretchedCoordinates(testCase)
            background = IMBackground.exponentialStratification(N0=.01,b=700,zDomain=[-1000 0],surfaceCondition="rigidLid");
            evp = IMInternalModes.waveModesAtWavenumber(background,k=1e-3);
            basis = IMSolverSpectral(nEVP=32,coordinateKind="wkb").solveEVP(evp,nModes=4);
            for variable = ["F","G"]
                normalized = basis.withNormalization("rms"+variable);
                measuredRms = sqrt(integral(@(z) abs(normalized.evaluate(z,variable=variable)).^2,-1000,0,ArrayValued=true,AbsTol=1e-10,RelTol=1e-12)/1000);
                testCase.verifyEqual(measuredRms,ones(1,4),AbsTol=3e-11);
            end
        end

        function canonicalAndRetiredNamesCannotBeRepurposed(testCase)
            bg = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            neutral = IMBackground(N2=@(z) zeros(size(z)),zDomain=bg.zDomain,surfaceCondition="freeSurface");
            critical = IMBackground(N2=@(z) bg.f0^2*ones(size(z)),zDomain=bg.zDomain,surfaceCondition="freeSurface");
            endpoint = IMGeostrophicZeroAPVModes.atWavenumber(bg,1e-3);
            bases = {IMSolverSpectral(nEVP=32).solveEVP(IMEigenvalueProblem(zDomain=[-1 0],p=1,q=0,r=1),nModes=2), ...
                createIMTestBasis("wave",bg,nModes=2),createIMTestBasis("wave",bg,nModes=2,solver=IMSolverAnalytical()), ...
                IMSolverSpectral(nEVP=32).solveGeostrophicZeroAPVModes(endpoint),IMSolverAnalytical().solveGeostrophicZeroAPVModes(endpoint), ...
                IMNeutralLabelModesBasis(neutral,neutralIntervals=[-1000 0],nModesPerInterval=2), ...
                IMZeroDepthWaveModesBasis(critical,1e-3,criticalIntervals=[-1000 0],nModesPerInterval=2)};
            canonical = ["rms","rmsF","rmsG","maxAbsF","maxAbsG","surfaceF","waveEnergy","generalizedPotentialEnstrophy","meanDensityAnomalyEnergy","endpointResponse"];
            retired = ["unity","kConstant","uMax","wMax","surfacePressure","geostrophic","depth","depthRMS","depthRMSF","depthRMSG","generalizedEnergy","innerProductG"];
            for j = 1:numel(bases)
                basis = bases{j}; prefix = string(class(basis));
                if isa(basis,"IMBasisSet"), prefix = "IMBasisSet"; end
                for name = canonical
                    testCase.verifyError(@() basis.addNormalization(name,@(~,~) 2),prefix+":ReservedNormalization");
                end
                for name = retired
                    testCase.verifyFalse(ismember(name,basis.normalizationNames()));
                    testCase.verifyError(@() basis.withNormalization(name),prefix+":UnsupportedNormalization");
                    testCase.verifyError(@() basis.addNormalization(name,@(~,~) 2),prefix+":UnsupportedNormalization");
                end
                testCase.verifyError(@() basis.withNormalization(Normalization.generalizedPotentialEnstrophy),prefix+":InvalidNormalizationName");
            end
        end

        function customReplacementDetachesCachesAndRetainsColumnIdentity(testCase)
            bg = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            for solver = {IMSolverSpectral(nEVP=40),IMSolverAnalytical()}
                original = createIMTestBasis("wave",bg,nModes=4,solver=solver{1});
                original.normalizationFactors();
                custom = original.addNormalization("customScale",@(~,j) j+1).withNormalization("customScale");
                custom.normalizationFactors();
                selected = custom.selectColumns([4 2]);
                replacement = custom.addNormalization("customScale",@(~,j) 2*(j+1));
                testCase.verifyEqual(selected.normalizationFactors(),[5 3]);
                testCase.verifyEqual(replacement.normalizationFactors(),2*(2:5));
                testCase.verifyEqual(custom.normalizationFactors(),2:5);
                values = custom.F([-800;-200]);
                testCase.verifyEqual(selected.F([-800;-200]),values(:,[4 2]),AbsTol=1e-12);
                testCase.verifyEqual(original.normalization,"waveEnergy");
                testCase.verifyEqual(selected.normSignatures,[1 1]);
                testCase.verifyEmpty(selected.selectColumns([]).normalizationFactors());
            end
        end

        function intervalCountsAreRequiredBeforeStratificationEvaluation(testCase)
            calls = 0;
            bg = IMBackground(N2=@sample,zDomain=[-1000 0],surfaceCondition="freeSurface");
            calls = 0;
            testCase.verifyError(@() IMNeutralLabelModesBasis(bg,neutralIntervals=[-1000 0]),"IMNeutralLabelModesBasis:MissingModeCount");
            testCase.verifyError(@() IMZeroDepthWaveModesBasis(bg,1e-3,criticalIntervals=[-1000 0]),"IMZeroDepthWaveModesBasis:MissingModeCount");
            testCase.verifyError(@() IMNeutralLabelModesBasis(bg,neutralIntervals=[-1000 0],nModesPerInterval=[]),"IMNeutralLabelModesBasis:MissingModeCount");
            testCase.verifyEqual(calls,0);
            basis = IMNeutralLabelModesBasis([-1000 -800;-200 0],zDomain=[-1000 0],nModesPerInterval=2);
            testCase.verifyEqual(basis.nModesPerInterval,[2 2]);
            testCase.verifyEqual(basis.selectColumns([4 2]).nModesPerInterval,[1 1]);
            testCase.verifyEqual(basis.constructionRecord.selection.request.nModesPerInterval,2);
            empty = IMNeutralLabelModesBasis([-1000 0],zDomain=[-1000 0],nModesPerInterval=0);
            testCase.verifyEqual(empty.nColumns,0);
            function values = sample(z)
                calls = calls+1; values = zeros(size(z));
            end
        end
    end
end
