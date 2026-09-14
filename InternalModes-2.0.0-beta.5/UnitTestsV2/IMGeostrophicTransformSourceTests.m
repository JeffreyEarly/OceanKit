classdef IMGeostrophicTransformSourceTests < matlab.unittest.TestCase

    properties (Access = private)
        originalPath
    end

    methods (TestClassSetup)
        function configurePath(testCase)
            repoRoot = fileparts(fileparts(mfilename("fullpath")));
            testCase.originalPath = path;
            addpath(repoRoot);
        end
    end

    methods (TestClassTeardown)
        function restorePath(testCase)
            path(testCase.originalPath);
        end
    end

    methods (Test)
        function constantSourceProjectionMatchesIndependentFunctionals(testCase)
            data = IMGeostrophicTransformSourceTests.buildCase();
            [vorticity,displacement] = IMGeostrophicTransformSourceTests.sourceFields(data.z,data.k,2);
            [Sq,S0] = data.transform.transformSourceForward(vorticitySource=vorticity,displacementSource=displacement);
            [expectedSq,expectedS0] = IMGeostrophicTransformSourceTests.independentProjection(data,vorticity,displacement);
            testCase.verifyEqual(Sq,expectedSq,RelTol=2e-12,AbsTol=2e-12)
            testCase.verifyEqual(S0,expectedS0,RelTol=3e-11,AbsTol=3e-11)
            testCase.verifyGreaterThan(norm(data.apvTransform.targetGramMatrix(variable="F")-eye(data.nModes),"fro"),1)
            testCase.verifyGreaterThan(norm(data.apvTransform.targetGramMatrix(variable="G")-eye(data.nModes),"fro"),1)
        end

        function rotatedSourceCoordinatesPreserveCanonicalProjection(testCase)
            data = IMGeostrophicTransformSourceTests.buildCase();
            [vorticity,displacement] = IMGeostrophicTransformSourceTests.sourceFields(data.z,data.k,3);
            rotated = data.zeroAPVModes.rotateBoundaryDepth(g0=data.g0,gd=data.gd);
            [SqCanonical,S0Canonical] = data.transform.transformSourceForward(vorticitySource=vorticity,displacementSource=displacement);
            [SqRotated,S0Rotated] = data.transform.transformSourceForward(vorticitySource=vorticity,displacementSource=displacement,zeroAPVCoordinates=rotated);
            testCase.verifyEqual(SqRotated,SqCanonical,AbsTol=0)
            for iK = 1:numel(data.k)
                expected = rotated.rotationMatrix(:,:,iK)\reshape(S0Canonical(:,iK,:),2,[]);
                testCase.verifyEqual(reshape(S0Rotated(:,iK,:),2,[]),expected,RelTol=3e-12,AbsTol=3e-12)
                testCase.verifyEqual(rotated.rotationMatrix(:,:,iK)*reshape(S0Rotated(:,iK,:),2,[]),reshape(S0Canonical(:,iK,:),2,[]),RelTol=3e-12,AbsTol=3e-12)
            end
        end

        function exponentialSourceProjectionMatchesIndependentFunctionals(testCase)
            N0 = 5.2e-3;
            b = 1300;
            N2 = @(z) N0^2*exp(2*z/b);
            data = IMGeostrophicTransformSourceTests.buildCase(N2=N2,zDomain=[-4000 0],k=[8e-6 1.5e-5],nEVP=160,nGrid=129);
            [vorticity,displacement] = IMGeostrophicTransformSourceTests.sourceFields(data.z,data.k,2);
            [Sq,S0] = data.transform.transformSourceForward(vorticitySource=vorticity,displacementSource=displacement);
            [expectedSq,expectedS0] = IMGeostrophicTransformSourceTests.independentProjection(data,vorticity,displacement);
            testCase.verifyEqual(Sq,expectedSq,RelTol=3e-11,AbsTol=3e-11)
            testCase.verifyEqual(S0,expectedS0,RelTol=2e-9,AbsTol=2e-9)
        end

        function sourceProjectionRequiresBothAPVChannels(testCase)
            data = IMGeostrophicTransformSourceTests.buildCase(variables="F");
            [vorticity,displacement] = IMGeostrophicTransformSourceTests.sourceFields(data.z,data.k,1);
            try
                data.transform.transformSourceForward(vorticitySource=vorticity,displacementSource=displacement);
                testCase.assertFail("Expected a missing G-channel failure.")
            catch exception
                testCase.verifyEqual(string(exception.identifier),"IMGeostrophicTransform:UnavailableAPVChannel")
                testCase.verifySubstring(exception.message,"not requested")
            end
        end

        function freeSurfaceSingularLimitRequiresSurfaceSourceSample(testCase)
            D = 1000;
            g = 9.81;
            N0 = 5.2e-3;
            N2 = @(z) N0^2*ones(size(z));
            g0 = -g;
            gd = Inf;
            k = 1e-4;
            solver = IMSolverSpectral(nEVP=112);
            evp = IMInternalModes.geostrophicAPVModes(N2=N2,zDomain=[-D 0],g=g,g0=g0,gd=gd);
            basis = solver.solveEVP(evp,nModes=8);
            z = linspace(-D,-D/64,64).';
            weights = ones(size(z))*(D/numel(z));
            apvTransform = basis.discreteTransform(z=z,weights=weights,nModes=5,variables=["F","G"],gramTolerance=100);
            problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=[-D 0],f0=1e-4,g=g,k=k,endpoints="surface");
            zeroAPVModes = solver.solveGeostrophicZeroAPVModes(problem);
            transform = IMGeostrophicTransform(apvTransform=apvTransform,zeroAPVModes=zeroAPVModes,g0=g0,gd=gd);
            testCase.verifyError(@() transform.transformSourceForward(vorticitySource=ones(numel(z),1),displacementSource=ones(numel(z),1)),"IMGeostrophicTransform:MissingSourceEndpointSample")
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
                options.variables string = ["F","G"]
                options.nModes (1,1) double = 6
                options.nEVP (1,1) double = 112
                options.nGrid (1,1) double = 65
            end
            if isempty(options.N2)
                N0 = 5.2e-3;
                N2 = @(z) N0^2*ones(size(z));
            else
                N2 = options.N2;
            end
            solver = IMSolverSpectral(nEVP=options.nEVP);
            evp = IMInternalModes.geostrophicAPVModes(N2=N2,zDomain=options.zDomain,g=options.g,g0=options.g0,gd=options.gd);
            basis = solver.solveEVP(evp,nModes=options.nModes+3);
            z = linspace(options.zDomain(1),options.zDomain(2),options.nGrid).';
            depth = diff(options.zDomain);
            weights = [0.5;ones(options.nGrid-2,1);0.5]*(depth/(options.nGrid-1));
            weightVariation = 1+0.12*sin(pi*(z-options.zDomain(1))/depth);
            weights = weights.*weightVariation;
            weights = weights*(depth/sum(weights));
            apvTransform = basis.discreteTransform(z=z,weights=weights,nModes=options.nModes,variables=options.variables,gramTolerance=100);
            problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=options.zDomain,f0=options.f0,g=options.g,k=options.k);
            zeroAPVModes = solver.solveGeostrophicZeroAPVModes(problem);
            transform = IMGeostrophicTransform(apvTransform=apvTransform,zeroAPVModes=zeroAPVModes,g0=options.g0,gd=options.gd);
            data = struct(N2=N2,zDomain=options.zDomain,g=options.g,f0=options.f0,g0=options.g0,gd=options.gd,k=reshape(options.k,1,[]), ...
                z=z,weights=weights,nModes=options.nModes,apvTransform=apvTransform,zeroAPVModes=zeroAPVModes,transform=transform);
        end

        function [vorticity,displacement] = sourceFields(z,k,nFields)
            nK = numel(k);
            depth = max(z)-min(z);
            vorticity = zeros(numel(z),nK,nFields);
            displacement = zeros(numel(z),nK,nFields);
            x = (z-min(z))/depth;
            for iField = 1:nFields
                for iK = 1:nK
                    phase = exp(1i*(0.2*iField+0.15*iK));
                    vorticity(:,iK,iField) = phase*(sin(pi*iField*x)+0.2*cos(2*pi*x));
                    displacement(:,iK,iField) = conj(phase)*(0.3+cos(pi*x)+0.1*iK*sin(3*pi*x));
                end
            end
        end

        function [Sq,S0] = independentProjection(data,vorticity,displacement)
            nK = numel(data.k);
            nFields = size(vorticity,3);
            nEndpoints = numel(data.zeroAPVModes.endpoints);
            depth = diff(data.zDomain);
            FInverse = data.apvTransform.inverseMatrix(variable="F");
            GInverse = data.apvTransform.inverseMatrix(variable="G");
            FMetric = data.apvTransform.metricMatrix(variable="F");
            GMetric = data.apvTransform.metricMatrix(variable="G");
            Sq = zeros(data.nModes,nK,nFields,"like",vorticity+displacement(1,1,1));
            S0 = zeros(nEndpoints,nK,nFields,"like",vorticity+displacement(1,1,1));
            zeroF = data.zeroAPVModes.F(data.z);
            zeroG = data.zeroAPVModes.G(data.z);
            endpointZ = [data.zDomain(2);data.zDomain(1)];
            zeroFEndpoints = data.zeroAPVModes.F(endpointZ);
            zeroGEndpoints = data.zeroAPVModes.G(endpointZ);
            generalizedEnergy = data.zeroAPVModes.generalizedEnergyMatrix(g0=data.g0,gd=data.gd);

            for iK = 1:nK
                omega = reshape(vorticity(:,iK,:),numel(data.z),nFields);
                eta = reshape(displacement(:,iK,:),numel(data.z),nFields);
                Sq(:,iK,:) = reshape((FInverse.'*FMetric*omega)/depth-(data.f0/depth)*(GInverse.'*GMetric*eta),data.nModes,1,nFields);

                FPairing = zeroF(:,:,iK).'*(data.weights.*omega)/depth;
                GPairing = zeroG(:,:,iK).'*((data.weights.*data.N2(data.z)/data.g).*eta);
                surfaceResponse = zeroGEndpoints(1,:,iK)-zeroFEndpoints(1,:,iK);
                GPairing = GPairing+(data.g0/data.g)*surfaceResponse.'*eta(end,:);
                GPairing = GPairing+(data.gd/data.g)*zeroGEndpoints(2,:,iK).'*eta(1,:);
                p0 = FPairing-(data.f0/depth)*GPairing;
                normalizedGram = (2*data.k(iK)^2/depth)*generalizedEnergy(:,:,iK);
                S0(:,iK,:) = reshape(normalizedGram\p0,nEndpoints,1,nFields);
            end
        end
    end
end
