classdef IMInternalModesDiscreteTransformTests < matlab.unittest.TestCase

    properties (Access = private)
        originalPath
        basisSet
        z
        weights
    end

    methods (TestClassSetup)
        function configureConstantFamily(testCase)
            repoRoot = fileparts(fileparts(mfilename("fullpath")));
            testCase.originalPath = path;
            addpath(repoRoot);
            D = 1000;
            N2 = @(z) 1e-4*ones(size(z));
            solver = IMSolverSpectral(nEVP=96);
            evp = IMInternalModes.hydrostaticFModes(N2=N2,zDomain=[-D 0]);
            testCase.basisSet = solver.solveEVP(evp,nModes=10);
            testCase.z = linspace(-D,0,21).';
            testCase.weights = [0.5;ones(19,1);0.5]*(D/20);
        end
    end

    methods (TestClassTeardown)
        function restorePath(testCase)
            path(testCase.originalPath);
        end
    end

    methods (Test)
        function alignedChannelsShareFamilyAndObeyActiveProjectors(testCase)
            variables = ["F","G"];
            transform = testCase.basisSet.discreteTransform(z=testCase.z,weights=testCase.weights,nModes=6,variables=variables,gramTolerance=1);

            testCase.verifyClass(transform,"IMInternalModesDiscreteTransform")
            testCase.verifyEqual(transform.availableVariables,variables)
            testCase.verifyEqual(transform.primaryVariable,"F")
            testCase.verifyEqual(transform.modeNumber,0:5,AbsTol=0)
            testCase.verifyEqual(transform.h,testCase.basisSet.h(1:6),AbsTol=0)
            expectedF = testCase.basisSet.F(testCase.z);
            expectedG = testCase.basisSet.G(testCase.z);
            testCase.verifyEqual(transform.inverseMatrix(variable="F"),expectedF(:,1:6),RelTol=1e-12,AbsTol=1e-12)
            testCase.verifyEqual(transform.inverseMatrix(variable="G"),expectedG(:,1:6),RelTol=1e-12,AbsTol=1e-12)
            testCase.verifyEqual(transform.activeModeMask(variable="F"),true(1,6))
            testCase.verifyEqual(transform.activeModeMask(variable="G"),[false true(1,5)])
            for variable = variables
                projector = diag(double(transform.activeModeMask(variable=variable)));
                testCase.verifyEqual(transform.forwardMatrix(variable=variable)*transform.inverseMatrix(variable=variable),projector,RelTol=1e-10,AbsTol=1e-11)
            end
            coefficients = reshape(1:12,6,2);
            expectedG = testCase.basisSet.G(testCase.z);
            testCase.verifyEqual(transform.transformBack(coefficients,variable="G"),expectedG(:,1:6)*coefficients,RelTol=1e-12,AbsTol=1e-12)
            testCase.verifyEqual(transform.transformBack(coefficients),transform.inverseMatrix(variable="F")*coefficients,RelTol=1e-12,AbsTol=1e-12)
        end

        function solvedGFormIsThePrimaryTransformVariable(testCase)
            evp = IMInternalModes.hydrostaticGModes(N2=@(z) 1e-4*ones(size(z)),zDomain=[-1000 0]);
            basis = IMSolverSpectral(nEVP=96).solveEVP(evp,nModes=6);
            transform = basis.discreteTransform(z=testCase.z,weights=testCase.weights,nModes=4,variables=["F","G"],gramTolerance=1);
            coefficients = reshape(1:8,4,2);

            testCase.verifyEqual(transform.primaryVariable,"G")
            testCase.verifyEqual(transform.inverseMatrix(),transform.inverseMatrix(variable="G"),AbsTol=0)
            testCase.verifyEqual(transform.transformBack(coefficients),transform.transformBack(coefficients,variable="G"),AbsTol=0)
        end

        function preparedPrefixAndReweightingMatchDirectConstruction(testCase)
            variables = ["F","G"];
            candidate = testCase.basisSet.discreteTransform(z=testCase.z,weights=testCase.weights,nModes=6,variables=variables,gramTolerance=1);

            prefix = candidate.prefixTransform(4);
            directPrefix = testCase.basisSet.discreteTransform(z=testCase.z,weights=testCase.weights,nModes=4,variables=variables,gramTolerance=1);
            testCase.verifyEqual(prefix.modeNumber,directPrefix.modeNumber,AbsTol=0)
            testCase.verifyEqual(prefix.h,directPrefix.h,AbsTol=0)
            for variable = variables
                testCase.verifyEqual(prefix.activeModeMask(variable=variable),directPrefix.activeModeMask(variable=variable))
                testCase.verifyEqual(prefix.inverseMatrix(variable=variable),directPrefix.inverseMatrix(variable=variable),RelTol=1e-13,AbsTol=1e-13)
                testCase.verifyEqual(prefix.endpointValues(variable=variable),directPrefix.endpointValues(variable=variable),RelTol=1e-13,AbsTol=1e-13)
                testCase.verifyEqual(prefix.metricMatrix(variable=variable),directPrefix.metricMatrix(variable=variable),RelTol=1e-13,AbsTol=1e-13)
                testCase.verifyEqual(prefix.targetGramMatrix(variable=variable),directPrefix.targetGramMatrix(variable=variable),RelTol=1e-12,AbsTol=1e-12)
                testCase.verifyEqual(prefix.forwardMatrix(variable=variable),directPrefix.forwardMatrix(variable=variable),RelTol=1e-11,AbsTol=1e-12)
            end

            newWeights = testCase.weights .* linspace(0.9,1.1,length(testCase.weights)).';
            reweighted = candidate.transformWithWeights(newWeights);
            directReweighted = testCase.basisSet.discreteTransform(z=testCase.z,weights=newWeights,nModes=6,variables=variables,gramTolerance=1);
            testCase.verifyEqual(reweighted.weights,newWeights,AbsTol=0)
            for variable = variables
                testCase.verifyEqual(reweighted.activeModeMask(variable=variable),directReweighted.activeModeMask(variable=variable))
                testCase.verifyEqual(reweighted.inverseMatrix(variable=variable),directReweighted.inverseMatrix(variable=variable),RelTol=1e-13,AbsTol=1e-13)
                testCase.verifyEqual(reweighted.metricMatrix(variable=variable),directReweighted.metricMatrix(variable=variable),RelTol=1e-13,AbsTol=1e-13)
                testCase.verifyEqual(reweighted.targetGramMatrix(variable=variable),directReweighted.targetGramMatrix(variable=variable),RelTol=1e-12,AbsTol=1e-12)
                testCase.verifyEqual(reweighted.gramMatrix(variable=variable),directReweighted.gramMatrix(variable=variable),RelTol=1e-12,AbsTol=1e-12)
                testCase.verifyEqual(reweighted.forwardMatrix(variable=variable),directReweighted.forwardMatrix(variable=variable),RelTol=1e-11,AbsTol=1e-12)
            end
        end

        function projectionFunctionalAndForwardSolveAreDistinct(testCase)
            transform = testCase.basisSet.discreteTransform(z=testCase.z,weights=testCase.weights,nModes=5,variables=["F","G"],gramTolerance=1);
            values = [sin(testCase.z/300) cos(testCase.z/210)];
            for variable = ["F","G"]
                inverse = transform.inverseMatrix(variable=variable);
                metric = transform.metricMatrix(variable=variable);
                active = transform.activeModeMask(variable=variable);
                expectedPairings = inverse.'*metric*values;
                expectedPairings(~active,:) = 0;
                testCase.verifyEqual(transform.modeProjectionFunctional(values,variable=variable),expectedPairings,RelTol=1e-12,AbsTol=1e-12)
                testCase.verifyEqual(transform.transformForward(values,variable=variable),transform.forwardMatrix(variable=variable)*values,RelTol=1e-12,AbsTol=1e-12)
            end
        end

        function endpointTracesAndPhysicalSnapshotAreDetached(testCase)
            transform = testCase.basisSet.discreteTransform(z=testCase.z,weights=testCase.weights,nModes=4,variables=["F","G"],gramTolerance=1);
            expectedF = [testCase.basisSet.F(0);testCase.basisSet.F(-1000)];
            expectedG = [testCase.basisSet.G(0);testCase.basisSet.G(-1000)];
            testCase.verifyEqual(transform.endpointLocations,["surface";"bottom"])
            testCase.verifyEqual(transform.endpointValues(variable="F"),expectedF(:,1:4),RelTol=1e-12,AbsTol=1e-12)
            testCase.verifyEqual(transform.endpointValues(variable="G"),expectedG(:,1:4),RelTol=1e-12,AbsTol=1e-12)
            testCase.verifyEqual(transform.N2Values,testCase.basisSet.N2(testCase.z),AbsTol=0)
            testCase.verifyEqual(transform.depth,1000,AbsTol=0)
            testCase.verifyEqual(transform.modeFamily,"hydrostatic")
            testCase.verifyFalse(isprop(transform,"basisSet"))
        end

        function unavailableContinuousAndCompanionEndpointChannelsExplainWhy(testCase)
            wave = IMInternalModes.waveModesAtFrequency(N2=@(z) 1e-4*ones(size(z)),zDomain=[-1000 0],omega=2e-3);
            waveSolver = IMSolverSpectral(nEVP=96);
            waveBasis = waveSolver.solveEVP(wave,nModes=3);
            testCase.verifyError(@() waveBasis.discreteTransform(z=testCase.z,weights=testCase.weights,nModes=2,variables="F"), ...
                "IMInternalModesBasis:UnavailableDiscreteTransformVariable")
            waveTransform = waveBasis.discreteTransform(z=testCase.z,weights=testCase.weights,nModes=2,variables="G",gramTolerance=1);
            testCase.verifyFalse(waveTransform.hasForwardTransform(variable="F"))
            testCase.verifySubstring(waveTransform.forwardTransformReason(variable="F"),"continuous inner product")
            testCase.verifySize(waveTransform.inverseMatrix(variable="F"),[length(testCase.z) 2])

            linked = IMInternalModes.hydrostaticFModes(N2=@(z) 1e-4*ones(size(z)),zDomain=[-1000 0], ...
                surfaceBoundary=IMBoundaryCondition(a=2,b=4,c=1,d=0));
            linkedSolver = IMSolverSpectral(nEVP=96);
            linkedBasis = linkedSolver.solveEVP(linked,nModes=3);
            try
                linkedBasis.discreteTransform(z=testCase.z,weights=testCase.weights,nModes=2,variables="G");
                testCase.assertFail("Expected a companion-variable endpoint error.")
            catch exception
                testCase.verifyEqual(string(exception.identifier),"IMInternalModesBasis:UnavailableDiscreteTransformVariable")
                testCase.verifySubstring(exception.message,"companion-variable")
            end
        end

        function apvMetadataPreservesEndpointParameters(testCase)
            g0 = -3.2;
            gd = 7.5;
            evp = IMInternalModes.geostrophicAPVModes(N2=@(z) 1e-4*ones(size(z)),zDomain=[-1000 0],g0=g0,gd=gd,surfaceBoundary="rigidLid");
            solver = IMSolverSpectral(nEVP=128);
            basis = solver.solveEVP(evp,nModes=4);
            transform = basis.discreteTransform(z=testCase.z,weights=testCase.weights,nModes=3,variables="F",gramTolerance=10);
            testCase.verifyEqual(transform.problemMetadata.g0,g0,AbsTol=0)
            testCase.verifyEqual(transform.problemMetadata.gd,gd,AbsTol=0)
            testCase.verifyEqual(transform.problemMetadata.surfaceBoundary,"rigidLid")
        end

        function sameVariableEndpointMetricsAreRepresentedAndDerivativeMetricsAreRejected(testCase)
            N2 = @(z) 1e-4*ones(size(z));
            g = 9.81;
            sameVariableEVP = IMInternalModes.hydrostaticFModes(N2=N2,zDomain=[-1000 0],g=g, ...
                surfaceBoundary=IMBoundaryCondition(a=2,b=3));
            solver = IMSolverSpectral(nEVP=96);
            basis = solver.solveEVP(sameVariableEVP,nModes=3);
            transform = basis.discreteTransform(z=testCase.z,weights=testCase.weights,nModes=2,variables="G",gramTolerance=10);
            expectedMetric = diag((N2(testCase.z)/g).*testCase.weights);
            expectedMetric(end,end) = expectedMetric(end,end)-3/(2*g);
            testCase.verifyEqual(transform.metricMatrix(variable="G"),expectedMetric,AbsTol=1e-13)

            derivativeEVP = IMInternalModes.hydrostaticFModes(N2=N2,zDomain=[-1000 0], ...
                surfaceBoundary=IMBoundaryCondition(a=1,b=0,c=0,d=1));
            derivativeBasis = solver.solveEVP(derivativeEVP,nModes=2);
            derivativeBasis.normalization = "unity";
            try
                derivativeBasis.discreteTransform(z=testCase.z,weights=testCase.weights,nModes=2,variables="F");
                testCase.assertFail("Expected a derivative endpoint representability error.")
            catch exception
                testCase.verifyEqual(string(exception.identifier),"IMInternalModesBasis:UnavailableDiscreteTransformVariable")
                testCase.verifySubstring(exception.message,"derivative trace")
            end
        end
    end
end
