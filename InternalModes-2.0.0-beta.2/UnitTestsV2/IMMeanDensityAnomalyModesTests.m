classdef IMMeanDensityAnomalyModesTests < matlab.unittest.TestCase

    properties (Access = private)
        originalPath
    end

    methods (TestClassSetup)
        function addRepositoryPath(testCase)
            repoRoot = fileparts(fileparts(mfilename("fullpath")));
            testCase.originalPath = path;
            addpath(repoRoot);
        end
    end

    methods (TestClassTeardown)
        function restoreRepositoryPath(testCase)
            path(testCase.originalPath);
        end
    end

    methods (Test)
        function factoryMapsEndpointLimitsAndStoresMetadata(testCase)
            [N2,zDomain,g] = testCase.constantProfile();
            cases = {
                0, 0, ["surface","bottom"]
                0.02, Inf, "surface"
                Inf, -0.03, "bottom"
                Inf, Inf, strings(1,0)};
            for iCase = 1:size(cases,1)
                g0 = cases{iCase,1};
                gd = cases{iCase,2};
                activeEndpoints = cases{iCase,3};
                evp = IMInternalModes.meanDensityAnomalyModes(N2=N2,zDomain=zDomain,g=g,g0=g0,gd=gd);

                testCase.verifyClass(evp,"IMMeanDensityAnomalyModes")
                testCase.verifyEqual(string(evp.name),"meanDensityAnomalyModes")
                testCase.verifyEqual(evp.formulation,"G")
                testCase.verifyEqual(evp.modeFamily,"meanDensityAnomaly")
                testCase.verifyEqual(evp.g0,g0,AbsTol=0)
                testCase.verifyEqual(evp.gd,gd,AbsTol=0)
                testCase.verifyEqual(evp.activeEndpoints,activeEndpoints)
                testCase.verifyEqual(evp.parameters.g0,g0,AbsTol=0)
                testCase.verifyEqual(evp.parameters.gd,gd,AbsTol=0)
                testCase.verifyEqual(evp.parameters.activeEndpoints,activeEndpoints)
                testCase.verifyEqual(testCase.boundaryCoefficients(evp.surfaceBoundary),testCase.surfaceCoefficients(g0,g),AbsTol=0)
                testCase.verifyEqual(testCase.boundaryCoefficients(evp.bottomBoundary),testCase.bottomCoefficients(gd,g),AbsTol=0)
            end
        end

        function endpointAccelerationsAreValidatedAndUnrelatedOptionsAreAbsent(testCase)
            [N2,zDomain] = testCase.constantProfile();
            testCase.verifyError(@() IMInternalModes.meanDensityAnomalyModes(N2=N2,zDomain=zDomain,g0=NaN,gd=1), ...
                "IMMeanDensityAnomalyModes:InvalidSurfaceAcceleration")
            testCase.verifyError(@() IMInternalModes.meanDensityAnomalyModes(N2=N2,zDomain=zDomain,g0=-Inf,gd=1), ...
                "IMMeanDensityAnomalyModes:InvalidSurfaceAcceleration")
            testCase.verifyError(@() IMInternalModes.meanDensityAnomalyModes(N2=N2,zDomain=zDomain,g0=1,gd=NaN), ...
                "IMMeanDensityAnomalyModes:InvalidBottomAcceleration")
            testCase.verifyError(@() IMInternalModes.meanDensityAnomalyModes(N2=N2,zDomain=zDomain,g0=1,gd=-Inf), ...
                "IMMeanDensityAnomalyModes:InvalidBottomAcceleration")
            testCase.verifyError(@() IMInternalModes.meanDensityAnomalyModes(N2=N2,zDomain=zDomain,g0=1,gd=1,f0=1e-4), ...
                "MATLAB:TooManyInputs")
            testCase.verifyError(@() IMInternalModes.meanDensityAnomalyModes(N2=N2,zDomain=zDomain,g0=1,gd=1,surfaceBoundary="freeSurface"), ...
                "MATLAB:TooManyInputs")
        end

        function generalizedEnergySpecificationIsExact(testCase)
            [N2,zDomain,g] = testCase.constantProfile();
            evp = IMInternalModes.meanDensityAnomalyModes(N2=N2,zDomain=zDomain,g=g,g0=-0.04,gd=0.03);
            GSpec = evp.innerProduct("G");
            FSpec = evp.innerProduct("F");
            z = linspace(zDomain(1),zDomain(2),7).';
            context = evp.contextForSolver(IMSolverSpectral(nEVP=32).configuredForEVP(evp));

            testCase.verifyTrue(GSpec.hasInnerProduct)
            testCase.verifyEqual(GSpec.interiorWeight(z,context),N2(z)/g,AbsTol=0)
            testCase.verifyEqual(GSpec.surfaceWeights.coefficient*GSpec.surfaceWeights.c^2,-0.04/g,RelTol=1e-14)
            testCase.verifyEqual(GSpec.bottomWeights.coefficient*GSpec.bottomWeights.c^2,0.03/g,RelTol=1e-14)
            testCase.verifyFalse(FSpec.hasInnerProduct)
            testCase.verifySubstring(FSpec.reason,"surface-referenced diagnostic pressure")
        end

        function constantModesAndIntegratedPressureAreExactInEverySpectralCoordinate(testCase)
            [N2,zDomain,g,N0] = testCase.constantProfile();
            D = diff(zDomain);
            z = linspace(zDomain(1),zDomain(2),401).';
            y = z-zDomain(1);
            for coordinateKind = ["z","wkb","density"]
                evp = IMInternalModes.meanDensityAnomalyModes(N2=N2,zDomain=zDomain,g=g,g0=0,gd=0);
                basis = IMSolverSpectral(nEVP=80,coordinateKind=coordinateKind).solveEVP(evp,nModes=5);
                G = basis.G(z);
                F = basis.F(z);
                for iMode = 1:5
                    n = basis.modeNumber(iMode);
                    amplitude = G(1,iMode);
                    if n == 0
                        expectedG = amplitude*ones(size(z));
                        expectedF = (N0^2/g)*amplitude*(zDomain(2)-z);
                    else
                        m = n*pi/D;
                        expectedG = amplitude*cos(m*y);
                        expectedF = -(N0^2/g)*(amplitude/m)*sin(m*y);
                    end
                    testCase.verifyEqual(G(:,iMode),expectedG,RelTol=2e-9,AbsTol=2e-10)
                    testCase.verifyEqual(F(:,iMode),expectedF,RelTol=2e-8,AbsTol=2e-10)
                end
                testCase.verifyEqual(F(end,:),zeros(1,5),AbsTol=2e-12)
            end
        end

        function endpointConfigurationsControlConstantNullMode(testCase)
            [N2,zDomain,g] = testCase.constantProfile();
            cases = {
                0.02, -0.01, true
                0.02, Inf, false
                Inf, -0.01, false
                Inf, Inf, false};
            for iCase = 1:size(cases,1)
                g0 = cases{iCase,1};
                gd = cases{iCase,2};
                expectNull = cases{iCase,3};
                evp = IMInternalModes.meanDensityAnomalyModes(N2=N2,zDomain=zDomain,g=g,g0=g0,gd=gd);
                basis = IMSolverSpectral(nEVP=96).solveEVP(evp,nModes=6);
                lambda = basis.eigenvalues;
                GSurface = basis.G(zDomain(2));
                GBottom = basis.G(zDomain(1));
                GzSurface = basis.Gz(zDomain(2));
                GzBottom = basis.Gz(zDomain(1));
                zCheck = linspace(zDomain(1),zDomain(2),501).';
                GScale = max(abs(basis.G(zCheck)),[],1);
                % A boundary-trapped mode can be tiny at the opposite end.
                % Include its full-domain derivative scale so roundoff in
                % that trace is not divided by a nearly zero local value.
                derivativeScale = GScale/diff(zDomain);

                testCase.verifyEqual(any(basis.modeNumber == 0),expectNull)
                if expectNull
                    zeroIndex = basis.modeNumber == 0;
                    testCase.verifyEqual(lambda(zeroIndex),0,AbsTol=0)
                    testCase.verifyEqual(basis.h(zeroIndex),Inf,AbsTol=0)
                end
                if isfinite(g0)
                    residual = GzSurface-(g0/g)*GSurface.*lambda;
                    residualScale = max([abs(GzSurface);abs((g0/g)*GSurface.*lambda);derivativeScale],[],1);
                    testCase.verifyLessThan(max(abs(residual)./max(residualScale,eps)),2e-8)
                else
                    testCase.verifyLessThan(max(abs(GSurface)./max(GScale,eps)),2e-9)
                end
                if isfinite(gd)
                    residual = GzBottom+(gd/g)*GBottom.*lambda;
                    residualScale = max([abs(GzBottom);abs((gd/g)*GBottom.*lambda);derivativeScale],[],1);
                    testCase.verifyLessThan(max(abs(residual)./max(residualScale,eps)),2e-8)
                else
                    testCase.verifyLessThan(max(abs(GBottom)./max(GScale,eps)),2e-9)
                end
            end
        end

        function signedNormalizationAndOrderingAreDeterministic(testCase)
            [N2,zDomain,g] = testCase.constantProfile();
            evp = IMInternalModes.meanDensityAnomalyModes(N2=N2,zDomain=zDomain,g=g,g0=-0.05,gd=0);
            solver = IMSolverSpectral(nEVP=96);
            first = solver.solveEVP(evp,nModes=6);
            second = solver.solveEVP(evp,nModes=6);
            gram = first.gramMatrix(variable="G");

            testCase.verifyEqual(first.modeNumber,[-1 0 1 2 3 4],AbsTol=0)
            testCase.verifyEqual(first.eigenvalues,second.eigenvalues,RelTol=2e-12,AbsTol=2e-12)
            testCase.verifyEqual(first.G(linspace(zDomain(1),zDomain(2),101).'), ...
                second.G(linspace(zDomain(1),zDomain(2),101).'),RelTol=2e-11,AbsTol=2e-11)
            testCase.verifyEqual(first.signatures,[-1 1 1 1 1 1],AbsTol=0)
            testCase.verifyEqual(gram,diag(first.signatures),RelTol=2e-8,AbsTol=2e-8)
            testCase.verifyEqual(string(first.normalization),"unity")
            testCase.verifyFalse(ismember("surfacePressure",first.normalizationNames()))
            testCase.verifyFalse(ismember("geostrophic",first.normalizationNames()))
        end

        function effectivelyInfinitePencilModeDoesNotDisplacePhysicalMode(testCase)
            [N2,zDomain,g] = testCase.constantProfile();
            evp = IMInternalModes.meanDensityAnomalyModes(N2=N2,zDomain=zDomain,g=g,g0=Inf,gd=-0.01);
            basis = IMSolverSpectral(nEVP=96).solveEVP(evp,nModes=4);

            testCase.verifyEqual(basis.modeNumber,[-1 1 2 3],AbsTol=0)
            testCase.verifyEqual(basis.eigenvalues(1),-g,RelTol=2e-7)
            testCase.verifyLessThan(max(abs(basis.G(zDomain(2)))./max(max(abs(basis.G(linspace(zDomain(1),zDomain(2),501).')),[],1),eps)),2e-9)
        end

        function exponentialPressureUsesSurfaceReferencedIntegration(testCase)
            D = 4000;
            g = 9.81;
            N0 = 5.2e-3;
            b = 1300;
            N2 = @(z) N0^2*exp(2*z/b);
            evp = IMInternalModes.meanDensityAnomalyModes(N2=N2,zDomain=[-D 0],g=g,g0=0.02,gd=Inf);
            basis = IMSolverSpectral(nEVP=128,coordinateKind="wkb").solveEVP(evp,nModes=8);
            z = linspace(-D,0,32001).';
            G = basis.G(z);
            integrand = (N2(z)/g).*G;
            expectedF = -flip(cumtrapz(flip(z),flip(integrand),1),1);
            F = basis.F(z);

            testCase.verifyLessThan(norm(F-expectedF,"fro")/norm(expectedF,"fro"),2e-7)
            testCase.verifyEqual(F(end,:),zeros(1,8),AbsTol=2e-11)
        end

        function finiteDifferencePressureUsesTheSameIntegralConvention(testCase)
            [N2,zDomain,g] = testCase.constantProfile();
            z = linspace(zDomain(1),zDomain(2),321).';
            evp = IMInternalModes.meanDensityAnomalyModes(N2=N2,zDomain=zDomain,g=g,g0=0,gd=0);
            basis = IMSolverFiniteDifference(z=z).solveEVP(evp,nModes=4);
            G = basis.G(z);
            expectedF = -flip(cumtrapz(flip(z),flip((N2(z)/g).*G),1),1);

            testCase.verifyEqual(basis.F(z),expectedF,RelTol=2e-11,AbsTol=2e-11)
            testCase.verifyEqual(basis.F(zDomain(2)),zeros(1,4),AbsTol=2e-12)
        end

        function singularConstantNullNormIsRejected(testCase)
            [N2,zDomain,g,N0] = testCase.constantProfile();
            singularG0 = -(N0^2*diff(zDomain));
            evp = IMInternalModes.meanDensityAnomalyModes(N2=N2,zDomain=zDomain,g=g,g0=singularG0,gd=0);
            testCase.verifyError(@() IMSolverSpectral(nEVP=64).solveEVP(evp,nModes=4), ...
                "IMMeanDensityAnomalyModesBasis:ZeroNormMode")
        end

        function discreteTransformProjectsGAndSynthesizesAlignedF(testCase)
            [N2,zDomain,g] = testCase.constantProfile();
            evp = IMInternalModes.meanDensityAnomalyModes(N2=N2,zDomain=zDomain,g=g,g0=0,gd=0);
            basis = IMSolverSpectral(nEVP=96).solveEVP(evp,nModes=18);
            [transform,assessment] = basis.discreteTransform(nPoints=16);
            coefficients = reshape(sin(1:(2*numel(transform.modeNumber))),numel(transform.modeNumber),2);
            GValues = transform.transformBack(coefficients,variable="G");

            testCase.verifyClass(transform,"IMInternalModesDiscreteTransform")
            testCase.verifyClass(assessment,"IMInternalModesDiscreteTransformAssessment")
            testCase.verifyEqual(transform.availableVariables,"G")
            testCase.verifyEqual(transform.primaryVariable,"G")
            testCase.verifyTrue(transform.hasForwardTransform(variable="G"))
            testCase.verifyFalse(transform.hasForwardTransform(variable="F"))
            testCase.verifyEqual(transform.transformForward(GValues,variable="G"),coefficients,RelTol=2e-10,AbsTol=2e-10)
            testCase.verifyEqual(transform.transformForward(GValues),coefficients,RelTol=2e-10,AbsTol=2e-10)
            testCase.verifyEqual(transform.transformBack(coefficients),GValues,RelTol=2e-10,AbsTol=2e-10)
            testCase.verifyEqual(transform.transformBack(coefficients,variable="F"), ...
                transform.inverseMatrix(variable="F")*coefficients,RelTol=2e-14,AbsTol=2e-14)
            testCase.verifyEqual(transform.endpointLocations,["surface";"bottom"])
            testCase.verifyEqual(transform.modeFamily,"meanDensityAnomaly")
            testCase.verifyEqual(transform.problemMetadata.g0,0,AbsTol=0)
            testCase.verifyEqual(transform.problemMetadata.gd,0,AbsTol=0)
            testCase.verifyEqual(transform.problemMetadata.activeEndpoints,["surface","bottom"])
            testCase.verifyError(@() transform.transformForward(transform.inverseMatrix(variable="F"),variable="F"), ...
                "IMInternalModesDiscreteTransform:UnavailableForwardTransform")
        end

        function sharedAPVGridRefitsMDAWeightsToTheSelectedBand(testCase)
            D = 4000;
            N0 = 5.2e-3;
            b = 1300;
            g = 9.81;
            Nz = 64;
            N2 = @(z) N0*N0*exp(2*z/b);
            solver = IMSolverSpectral(nEVP=256,coordinateKind="wkb");
            apv = solver.solveEVP(IMInternalModes.geostrophicAPVModes( ...
                N2=N2,zDomain=[-D 0],g=g,g0=0.02,gd=0.03,surfaceBoundary="freeSurface"),nModes=Nz+8);
            mda = solver.solveEVP(IMInternalModes.meanDensityAnomalyModes( ...
                N2=N2,zDomain=[-D 0],g=g,g0=0.02,gd=0.03),nModes=Nz+8);
            [z,apvGridDesign] = apv.modeRootGrid(nPoints=Nz);

            [transform,assessment] = mda.certifiedDiscreteTransform( ...
                z=z,variables="G",gridDesign=apvGridDesign,gramTolerance=1e-2);
            [~,allBandFit] = mda.quadratureWeightsForPoints(z=z,nModes=Nz,variables="G");
            [~,fixedRuleAssessment] = mda.discreteTransform( ...
                z=z,weights=allBandFit.transform.weights,variables="G",gramTolerance=1e-2);

            testCase.verifyGreaterThan(assessment.retainedModeCount,fixedRuleAssessment.retainedModeCount)
            testCase.verifyGreaterThanOrEqual(assessment.retainedModeCount,35)
            testCase.verifyEqual(assessment.weightFitModeCount,assessment.retainedModeCount)
            testCase.verifyEqual(length(transform.modeNumber),assessment.retainedModeCount)
            testCase.verifyEqual(assessment.gridDesign.sourceEVP,"geostrophicAPVModes")
            testCase.verifyEqual(assessment.gridDesign.sourceFamily,"hydrostatic")
            testCase.verifyEqual(assessment.gridDesign.generatingVariable,"F")
            testCase.verifyLessThanOrEqual(assessment.gramPolicy.error(end),1e-2)
            testCase.verifyTrue(assessment.certificationSearch.accepted(end))
        end

        function retainedPoliciesRespectSignedTargetsAndUndefinedProducts(testCase)
            [N2,zDomain,g] = testCase.constantProfile();
            positiveEVP = IMInternalModes.meanDensityAnomalyModes(N2=N2,zDomain=zDomain,g=g,g0=0,gd=0);
            positiveBasis = IMSolverSpectral(nEVP=96).solveEVP(positiveEVP,nModes=12);
            testCase.verifyError(@() positiveBasis.discreteTransform(nPoints=10,quadraticAliasingTolerance=0.1), ...
                "IMMeanDensityAnomalyModesBasis:UnavailableQuadraticAliasingPolicy")

            signedEVP = IMInternalModes.meanDensityAnomalyModes(N2=N2,zDomain=zDomain,g=g,g0=-0.05,gd=0);
            signedBasis = IMSolverSpectral(nEVP=96).solveEVP(signedEVP,nModes=12);
            testCase.verifyError(@() signedBasis.discreteTransform(nPoints=10,leakageTolerance=0.1,nCheckModes=11), ...
                "IMBasisSet:UnavailableDiscreteTransformPolicy")
        end
    end

    methods (Static, Access = private)
        function [N2,zDomain,g,N0] = constantProfile()
            N0 = 1e-2;
            zDomain = [-1000 0];
            g = 9.81;
            N2 = @(z) N0^2*ones(size(z));
        end

        function coefficients = boundaryCoefficients(boundary)
            coefficients = [boundary.a boundary.b boundary.c boundary.d];
        end

        function coefficients = surfaceCoefficients(g0,g)
            if isinf(g0)
                coefficients = [1 0 0 0];
            else
                coefficients = [0 1 g0/g 0];
            end
        end

        function coefficients = bottomCoefficients(gd,g)
            if isinf(gd)
                coefficients = [1 0 0 0];
            else
                coefficients = [0 1 -gd/g 0];
            end
        end
    end
end
