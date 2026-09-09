classdef IMGeostrophicTransformTests < matlab.unittest.TestCase

    properties (Access = private)
        originalPath
        standard
    end

    methods (TestClassSetup)
        function configureStandardCase(testCase)
            repoRoot = fileparts(fileparts(mfilename("fullpath")));
            testCase.originalPath = path;
            addpath(repoRoot);
            testCase.standard = IMGeostrophicTransformTests.buildCase();
        end
    end

    methods (TestClassTeardown)
        function restorePath(testCase)
            path(testCase.originalPath);
        end
    end

    methods (Test)
        function constructionCapturesCompatibilityAndEndpointResponse(testCase)
            data = testCase.standard;
            transform = data.transform;
            testCase.verifyEqual(transform.k,data.k,AbsTol=0)
            testCase.verifyEqual(transform.activeEndpoints,["surface","bottom"])
            testCase.verifyEqual(transform.g0,data.g0,AbsTol=0)
            testCase.verifyEqual(transform.gd,data.gd,AbsTol=0)
            testCase.verifyEqual(transform.f0,data.f0,AbsTol=0)
            testCase.verifySize(transform.apvEndpointResponse,[2 data.nModes numel(data.k)])

            FEndpoints = data.apvTransform.endpointValues(variable="F");
            GEndpoints = data.apvTransform.endpointValues(variable="G");
            boundaryValues = [GEndpoints(1,:)-FEndpoints(1,:);GEndpoints(2,:)];
            mu = data.k.^2+(data.f0^2/data.g)./data.apvTransform.h.';
            expected = zeros(size(transform.apvEndpointResponse));
            for iK = 1:numel(data.k)
                expected(:,:,iK) = boundaryValues.*reshape((-data.f0/data.g)./mu(:,iK),1,[]);
            end
            testCase.verifyEqual(transform.apvEndpointResponse,expected,RelTol=1e-13,AbsTol=1e-13)
            testCase.verifyEqual(transform.compatibilityDiagnostics.mu,mu,RelTol=1e-14)
            testCase.verifyEqual(transform.compatibilityDiagnostics.stratificationRelativeMismatch,0,AbsTol=0)
            testCase.verifyTrue(transform.compatibilityDiagnostics.endpointParameterMatches.g0)
            testCase.verifyTrue(transform.compatibilityDiagnostics.endpointParameterMatches.gd)
            testCase.verifyTrue(transform.compatibilityDiagnostics.surfaceConventionMatches)
        end

        function stateRoundTripPreservesComplexPagesAndTrailingDimensions(testCase)
            data = testCase.standard;
            Aq = randn(data.nModes,numel(data.k),2,3)+1i*randn(data.nModes,numel(data.k),2,3);
            A0 = randn(2,numel(data.k),2,3)+1i*randn(2,numel(data.k),2,3);
            [q,b] = data.transform.transformStateBack(APVCoefficients=Aq,zeroAPVCoefficients=A0);
            [AqBack,A0Back] = data.transform.transformStateForward(APV=q,endpointAnomalies=b);
            testCase.verifySize(q,[numel(data.z) numel(data.k) 2 3])
            testCase.verifySize(b,[2 numel(data.k) 2 3])
            testCase.verifyEqual(AqBack,Aq,RelTol=3e-11,AbsTol=3e-11)
            testCase.verifyEqual(A0Back,A0,RelTol=3e-10,AbsTol=3e-10)
        end

        function rotatedCoordinatesPreservePhysicalState(testCase)
            data = testCase.standard;
            rotated = data.zeroAPVModes.rotateBoundaryDepth(g0=data.g0,gd=data.gd);
            Aq = randn(data.nModes,numel(data.k),2)+1i*randn(data.nModes,numel(data.k),2);
            A0 = randn(2,numel(data.k),2)+1i*randn(2,numel(data.k),2);
            [qCanonical,bCanonical] = data.transform.transformStateBack(APVCoefficients=Aq,zeroAPVCoefficients=A0);
            [~,A0Rotated] = data.transform.transformStateForward(APV=qCanonical,endpointAnomalies=bCanonical,zeroAPVCoordinates=rotated);
            [qRotated,bRotated] = data.transform.transformStateBack(APVCoefficients=Aq,zeroAPVCoefficients=A0Rotated,zeroAPVCoordinates=rotated);
            testCase.verifyEqual(qRotated,qCanonical,RelTol=1e-13,AbsTol=1e-13)
            testCase.verifyEqual(bRotated,bCanonical,RelTol=3e-11,AbsTol=3e-11)
            for iK = 1:numel(data.k)
                expected = rotated.rotationMatrix(:,:,iK)\reshape(A0(:,iK,:),2,[]);
                testCase.verifyEqual(reshape(A0Rotated(:,iK,:),2,[]),expected,RelTol=2e-12,AbsTol=2e-12)
            end
        end

        function activeEndpointSubsetsIncludeFiniteZero(testCase)
            surface = IMGeostrophicTransformTests.buildCase(g0=0,gd=Inf,k=1.3e-4);
            bottom = IMGeostrophicTransformTests.buildCase(g0=Inf,gd=0,k=1.3e-4);
            testCase.verifyEqual(surface.transform.activeEndpoints,"surface")
            testCase.verifyEqual(bottom.transform.activeEndpoints,"bottom")

            AqSurface = randn(surface.nModes,1);
            A0Surface = randn(1,1);
            [qSurface,bSurface] = surface.transform.transformStateBack(APVCoefficients=AqSurface,zeroAPVCoefficients=A0Surface);
            [AqSurfaceBack,A0SurfaceBack] = surface.transform.transformStateForward(APV=qSurface,endpointAnomalies=bSurface);
            testCase.verifyEqual(AqSurfaceBack,AqSurface,RelTol=2e-11,AbsTol=2e-11)
            testCase.verifyEqual(A0SurfaceBack,A0Surface,RelTol=2e-11,AbsTol=2e-11)

            AqBottom = randn(bottom.nModes,1);
            A0Bottom = randn(1,1);
            [qBottom,bBottom] = bottom.transform.transformStateBack(APVCoefficients=AqBottom,zeroAPVCoefficients=A0Bottom);
            [AqBottomBack,A0BottomBack] = bottom.transform.transformStateForward(APV=qBottom,endpointAnomalies=bBottom);
            testCase.verifyEqual(AqBottomBack,AqBottom,RelTol=2e-11,AbsTol=2e-11)
            testCase.verifyEqual(A0BottomBack,A0Bottom,RelTol=2e-11,AbsTol=2e-11)
        end

        function rigidLidUsesGSurfaceResponse(testCase)
            data = IMGeostrophicTransformTests.buildCase(surfaceBoundary="rigidLid",k=1.2e-4);
            GEndpoints = data.apvTransform.endpointValues(variable="G");
            mu = data.k.^2+(data.f0^2/data.g)./data.apvTransform.h.';
            expectedSurface = GEndpoints(1,:).*((-data.f0/data.g)./mu(:,1).');
            testCase.verifyEqual(data.transform.apvEndpointResponse(1,:,1),expectedSurface,RelTol=1e-13,AbsTol=1e-13)
        end

        function constructorRejectsIncompatibleInputs(testCase)
            data = testCase.standard;
            testCase.verifyError(@() IMGeostrophicTransform(apvTransform=data.apvTransform,zeroAPVModes=data.zeroAPVModes,g0=data.g0+1e-3,gd=data.gd),"IMGeostrophicTransform:EndpointParameterMismatch")
            testCase.verifyError(@() IMGeostrophicTransform(apvTransform=data.apvTransform,zeroAPVModes=data.zeroAPVModes,g0=Inf,gd=Inf),"IMGeostrophicTransform:NoActiveEndpoint")
            rotated = data.zeroAPVModes.rotateBoundaryDepth(g0=data.g0,gd=data.gd);
            testCase.verifyError(@() IMGeostrophicTransform(apvTransform=data.apvTransform,zeroAPVModes=rotated,g0=data.g0,gd=data.gd),"IMGeostrophicTransform:CanonicalZeroAPVModesRequired")

            wrongEndpointProblem = IMGeostrophicZeroAPVModes.atWavenumber(N2=data.N2,zDomain=data.zDomain,f0=data.f0,g=data.g,k=data.k,endpoints="surface");
            wrongEndpoints = data.solver.solveGeostrophicZeroAPVModes(wrongEndpointProblem);
            testCase.verifyError(@() IMGeostrophicTransform(apvTransform=data.apvTransform,zeroAPVModes=wrongEndpoints,g0=data.g0,gd=data.gd),"IMGeostrophicTransform:ActiveEndpointMismatch")

            rigidProblem = IMGeostrophicZeroAPVModes.atWavenumber(N2=data.N2,zDomain=data.zDomain,f0=data.f0,g=data.g,k=data.k,surfaceBoundary="rigidLid");
            rigidModes = data.solver.solveGeostrophicZeroAPVModes(rigidProblem);
            testCase.verifyError(@() IMGeostrophicTransform(apvTransform=data.apvTransform,zeroAPVModes=rigidModes,g0=data.g0,gd=data.gd),"IMGeostrophicTransform:SurfaceConventionMismatch")

            wrongN2 = @(z) 1.01*data.N2(z);
            wrongN2Problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=wrongN2,zDomain=data.zDomain,f0=data.f0,g=data.g,k=data.k);
            wrongN2Modes = data.solver.solveGeostrophicZeroAPVModes(wrongN2Problem);
            testCase.verifyError(@() IMGeostrophicTransform(apvTransform=data.apvTransform,zeroAPVModes=wrongN2Modes,g0=data.g0,gd=data.gd),"IMGeostrophicTransform:StratificationMismatch")

            wrongDomainProblem = IMGeostrophicZeroAPVModes.atWavenumber(N2=data.N2,zDomain=[data.zDomain(1)+10 data.zDomain(2)],f0=data.f0,g=data.g,k=data.k);
            wrongDomainModes = data.solver.solveGeostrophicZeroAPVModes(wrongDomainProblem);
            testCase.verifyError(@() IMGeostrophicTransform(apvTransform=data.apvTransform,zeroAPVModes=wrongDomainModes,g0=data.g0,gd=data.gd),"IMGeostrophicTransform:DomainMismatch")

            wrongGravityProblem = IMGeostrophicZeroAPVModes.atWavenumber(N2=data.N2,zDomain=data.zDomain,f0=data.f0,g=data.g+0.1,k=data.k);
            wrongGravityModes = data.solver.solveGeostrophicZeroAPVModes(wrongGravityProblem);
            testCase.verifyError(@() IMGeostrophicTransform(apvTransform=data.apvTransform,zeroAPVModes=wrongGravityModes,g0=data.g0,gd=data.gd),"IMGeostrophicTransform:GravityMismatch")
            testCase.verifyError(@() IMGeostrophicTransform(apvTransform=data.apvTransform,zeroAPVModes=data.zeroAPVModes,g0=NaN,gd=data.gd),"IMGeostrophicTransform:InvalidEndpointAcceleration")
            testCase.verifyError(@() IMGeostrophicTransform(apvTransform=data.apvTransform,zeroAPVModes=data.zeroAPVModes,g0=data.g0,gd=-Inf),"IMGeostrophicTransform:InvalidEndpointAcceleration")
        end

        function methodValidationReportsChannelsCoordinatesAndShapes(testCase)
            data = testCase.standard;
            otherSolver = IMSolverSpectral(nEVP=80);
            otherCoordinates = otherSolver.solveGeostrophicZeroAPVModes(data.zeroAPVModes.problem);
            Aq = ones(data.nModes,numel(data.k));
            A0 = ones(2,numel(data.k));
            testCase.verifyError(@() data.transform.transformStateBack(APVCoefficients=Aq,zeroAPVCoefficients=A0,zeroAPVCoordinates=otherCoordinates),"IMGeostrophicTransform:IncompatibleZeroAPVCoordinates")
            testCase.verifyError(@() data.transform.transformStateBack(APVCoefficients=ones(data.nModes+1,numel(data.k)),zeroAPVCoefficients=A0),"IMGeostrophicTransform:InvalidArrayShape")
            testCase.verifyError(@() data.transform.transformStateForward(APV=ones(numel(data.z),numel(data.k),2),endpointAnomalies=ones(2,numel(data.k),3)),"IMGeostrophicTransform:TrailingDimensionMismatch")

            unavailable = IMGeostrophicTransformTests.buildCase(variables="G");
            q = ones(numel(unavailable.z),numel(unavailable.k));
            b = ones(2,numel(unavailable.k));
            try
                unavailable.transform.transformStateForward(APV=q,endpointAnomalies=b);
                testCase.assertFail("Expected the unavailable F channel to reject state projection.")
            catch exception
                testCase.verifyEqual(string(exception.identifier),"IMGeostrophicTransform:UnavailableAPVChannel")
                testCase.verifySubstring(exception.message,"not requested")
            end
        end

        function nearSingularMuIsRejected(testCase)
            D = 1000;
            g = 9.81;
            f0 = 1e-4;
            g0 = -0.01;
            gd = Inf;
            N2 = @(z) (5.2e-3)^2*ones(size(z));
            solver = IMSolverSpectral(nEVP=96);
            evp = IMInternalModes.geostrophicAPVModes(N2=N2,zDomain=[-D 0],g=g,g0=g0,gd=gd);
            basis = solver.solveEVP(evp,nModes=8);
            z = linspace(-D,0,65).';
            weights = [0.5;ones(63,1);0.5]*(D/64);
            apvTransform = basis.discreteTransform(z=z,weights=weights,nModes=6,variables=["F","G"],gramTolerance=100);
            negativeIndex = find(apvTransform.h < 0,1);
            testCase.assertNotEmpty(negativeIndex)
            singularK = sqrt(-f0^2/(g*apvTransform.h(negativeIndex)));
            problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=[-D 0],f0=f0,g=g,k=singularK,endpoints="surface");
            zeroAPVModes = solver.solveGeostrophicZeroAPVModes(problem);
            testCase.verifyError(@() IMGeostrophicTransform(apvTransform=apvTransform,zeroAPVModes=zeroAPVModes,g0=g0,gd=gd),"IMGeostrophicTransform:NearSingularMu")
        end

        function singularZeroAPVGeneralizedEnergyIsRejected(testCase)
            D = 1000;
            g = 9.81;
            f0 = 1e-4;
            gd = Inf;
            k = 1e-4;
            N2 = @(z) (5.2e-3)^2*ones(size(z));
            solver = IMSolverSpectral(nEVP=96);
            problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=[-D 0],f0=f0,g=g,k=k,endpoints="surface");
            zeroAPVModes = solver.solveGeostrophicZeroAPVModes(problem);
            g0 = -zeroAPVModes.energyMatrix/zeroAPVModes.surfaceBuoyancyMatrix;
            evp = IMInternalModes.geostrophicAPVModes(N2=N2,zDomain=[-D 0],g=g,g0=g0,gd=gd);
            basis = solver.solveEVP(evp,nModes=8);
            z = linspace(-D,0,65).';
            weights = [0.5;ones(63,1);0.5]*(D/64);
            apvTransform = basis.discreteTransform(z=z,weights=weights,nModes=6,variables=["F","G"],gramTolerance=100);
            channelData = struct();
            for variable = ["F","G"]
                channelData.(variable) = struct(available=true,reason="",activeModeMask=apvTransform.activeModeMask(variable=variable), ...
                    metricMatrix=apvTransform.metricMatrix(variable=variable),targetGramMatrix=apvTransform.targetGramMatrix(variable=variable), ...
                    targetMajorantGramMatrix=apvTransform.targetMajorantGramMatrix(variable=variable));
            end
            nonsingularAPV = IMInternalModesDiscreteTransform(z=apvTransform.z,weights=apvTransform.weights,modeNumber=apvTransform.modeNumber, ...
                h=abs(apvTransform.h)+1,normalization=apvTransform.normalization,inverseF=apvTransform.inverseMatrix(variable="F"), ...
                inverseG=apvTransform.inverseMatrix(variable="G"),endpointF=apvTransform.endpointValues(variable="F"), ...
                endpointG=apvTransform.endpointValues(variable="G"),channelData=channelData,zDomain=apvTransform.zDomain,g=apvTransform.g, ...
                modeFamily=apvTransform.modeFamily,N2Values=apvTransform.N2Values,problemMetadata=apvTransform.problemMetadata);
            testCase.verifyError(@() IMGeostrophicTransform(apvTransform=nonsingularAPV,zeroAPVModes=zeroAPVModes,g0=g0,gd=gd),"IMGeostrophicTransform:SingularZeroAPVGram")
        end

        function zeroWavenumberIsRejectedByProblemBoundary(testCase)
            N2 = @(z) ones(size(z));
            testCase.verifyError(@() IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=[-1 0],f0=1,k=0),"MATLAB:validators:mustBePositive")
        end
    end

    methods (Static, Access = private)
        function data = buildCase(options)
            arguments
                options.N2 = []
                options.zDomain (1,2) double = [-1000 0]
                options.g (1,1) double = 9.81
                options.f0 (1,1) double = 1e-4
                options.g0 (1,1) double = 0.03
                options.gd (1,1) double = 0.01
                options.k double = [1e-4 2e-4]
                options.surfaceBoundary string = "freeSurface"
                options.variables string = ["F","G"]
                options.nModes (1,1) double = 6
            end
            if isempty(options.N2)
                N0 = 5.2e-3;
                N2 = @(z) N0^2*ones(size(z));
            else
                N2 = options.N2;
            end
            solver = IMSolverSpectral(nEVP=96);
            evp = IMInternalModes.geostrophicAPVModes(N2=N2,zDomain=options.zDomain,g=options.g,g0=options.g0,gd=options.gd,surfaceBoundary=options.surfaceBoundary);
            basis = solver.solveEVP(evp,nModes=options.nModes+2);
            z = linspace(options.zDomain(1),options.zDomain(2),65).';
            depth = diff(options.zDomain);
            weights = [0.5;ones(63,1);0.5]*(depth/64);
            apvTransform = basis.discreteTransform(z=z,weights=weights,nModes=options.nModes,variables=options.variables,gramTolerance=10);
            canonicalEndpoints = ["surface","bottom"];
            endpoints = canonicalEndpoints([isfinite(options.g0),isfinite(options.gd)]);
            problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=options.zDomain,f0=options.f0,g=options.g,k=options.k,endpoints=endpoints,surfaceBoundary=options.surfaceBoundary);
            zeroAPVModes = solver.solveGeostrophicZeroAPVModes(problem);
            transform = IMGeostrophicTransform(apvTransform=apvTransform,zeroAPVModes=zeroAPVModes,g0=options.g0,gd=options.gd);
            data = struct(N2=N2,zDomain=options.zDomain,g=options.g,f0=options.f0,g0=options.g0,gd=options.gd,k=reshape(options.k,1,[]), ...
                solver=solver,z=z,weights=weights,nModes=options.nModes,apvTransform=apvTransform,zeroAPVModes=zeroAPVModes,transform=transform);
        end
    end
end
