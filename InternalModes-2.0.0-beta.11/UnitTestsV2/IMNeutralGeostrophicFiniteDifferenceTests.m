classdef IMNeutralGeostrophicFiniteDifferenceTests < matlab.unittest.TestCase
    properties (Access = private)
        originalPath
    end

    methods (TestClassSetup)
        function addRepositoryPath(testCase)
            repoRoot = fileparts(fileparts(mfilename("fullpath")));
            testCase.originalPath = path;
            addpath(repoRoot);
            addpath(fullfile(repoRoot,"UnitTestsV2","TestSupport"));
        end
    end

    methods (TestClassTeardown)
        function restorePath(testCase)
            path(testCase.originalPath);
        end
    end

    methods (Test)
        function constantAPVModesHaveNoCheckerboardBranchesAndConverge(testCase)
            N0 = 5e-3;
            evp = IMInternalModes.geostrophicAPVModes(N2=@(z) N0^2*ones(size(z)),zDomain=[-1000 0],g0=0,gd=0,surfaceBoundary="freeSurface");
            exact = N0^2*1000^2./(evp.g*((1:6)*pi).^2);
            errors = zeros(3,6);
            resolutions = [65 129 257];
            for i = 1:numel(resolutions)
                n = resolutions(i);
                solver = IMSolverFiniteDifference(z=linspace(-1000,0,n).');
                basis = solver.solveEVP(evp,nModes=6);
                testCase.verifyEqual(basis.modeNumber,1:6)
                testCase.verifySize(basis.nativeModes,[n 6])
                discrete = N0^2*(1000/(n-1))^2./(4*evp.g*sin((1:6)*pi/(2*(n-1))).^2);
                testCase.verifyEqual(basis.h,discrete,RelTol=2e-9)
                errors(i,:) = abs(basis.h./exact-1);
            end
            testCase.verifyLessThan(max(errors(1,:)),0.01)
            testCase.verifyGreaterThan(min(errors(1,:)./errors(2,:)),3.8)
            testCase.verifyGreaterThan(min(errors(2,:)./errors(3,:)),3.8)
        end

        function neutralCapFirstModesConvergeWithoutAliases(testCase)
            evp = IMInternalModes.geostrophicAPVModes(N2=@IMNeutralGeostrophicFiniteDifferenceTests.neutralCap,zDomain=[-1000 0],g0=0,gd=0,surfaceBoundary="freeSurface");
            reference = IMSolverSpectral(nEVP=512).solveEVP(evp,nModes=6);
            errors = zeros(3,6);
            for i = 1:3
                n = 2^(i+5)+1;
                basis = IMSolverFiniteDifference(z=linspace(-1000,0,n).').solveEVP(evp,nModes=6);
                testCase.verifyEqual(basis.modeNumber,1:6)
                errors(i,:) = abs(basis.h./reference.h-1);
                capF = basis.F(linspace(-180,-10,31).');
                testCase.verifyEqual(capF,zeros(size(capF)),AbsTol=1e-9)
                testCase.verifyGreaterThan(min(basis.inverseDepth),0)
            end
            testCase.verifyLessThan(max(errors(1,:)),0.04)
            testCase.verifyGreaterThan(min(errors(1,:)./errors(2,:)),3)
            testCase.verifyGreaterThan(min(errors(2,:)./errors(3,:)),3)
        end

        function signedNeutralCapSectorsMatchTheRegularSpectralReference(testCase)
            endpoints = [0,-0.01;-0.0175,0.0175];
            for j = 1:size(endpoints,1)
                evp = IMInternalModes.geostrophicAPVModes(N2=@IMNeutralGeostrophicFiniteDifferenceTests.neutralCap, ...
                    zDomain=[-1000 0],g0=endpoints(j,1),gd=endpoints(j,2),surfaceBoundary="freeSurface");
                reference = IMSolverSpectral(nEVP=512).solveEVP(evp,nModes=6);
                basis = IMSolverFiniteDifference(z=linspace(-1000,0,257).').solveEVP(evp,nModes=6);
                testCase.verifyEqual(basis.modeNumber,reference.modeNumber)
                testCase.verifyEqual(sign(basis.inverseDepth),sign(reference.inverseDepth))
                testCase.verifyLessThan(max(abs(basis.h./reference.h-1)),0.01)
            end
        end

        function enstrophyPublicPencilAndMimeticRayleighIdentityAgree(testCase)
            N2 = @IMNeutralGeostrophicFiniteDifferenceTests.neutralCap;
            evp = IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes(N2=N2,zDomain=[-1000 0],k=2*pi/1e4,latitude=45,alpha0=2e-10,alphaD=3e-10);
            solver = IMSolverFiniteDifference(z=linspace(-1000,0,129).');
            basis = solver.solveEVP(evp,nModes=6);
            z = solver.zNative;
            h = -diff(z);
            dual = [z(1);(z(1:end-1)+z(2:end))/2;z(end)];
            WF = [h(1)/2;(h(1:end-1)+h(2:end))/2;h(end)/2];
            F = basis.rawVariable("F",z);
            G = basis.rawVariable("G",dual);
            cellN2 = zeros(numel(h),1);
            for j = 1:numel(h), cellN2(j) = integral(N2,z(j+1),z(j),AbsTol=0,RelTol=1e-12); end
            volumeF = sum(WF.*F.^2,1);
            volumeG = sum(cellN2.*G(2:end-1,:).^2,1);
            numerator = evp.parameters.k^2*volumeF+evp.f0^2/evp.g^2*volumeG+evp.f0^2/evp.g*F(1,:).^2;
            denominator = volumeF+evp.f0^2/evp.parameters.alpha0*F(1,:).^2+evp.f0^2/evp.parameters.alphaD*F(end,:).^2;
            testCase.verifyEqual(basis.eigenvalues,numerator./denominator,RelTol=2e-12)
            testCase.verifyGreaterThan(min(basis.eigenvalues),0)
            [A,B] = evp.assemble(solver);
            V = [F;G];
            residual = vecnorm(A*V-(B*V).*basis.eigenvalues,2,1)./(norm(A,1)*vecnorm(V,2,1));
            testCase.verifyLessThan(max(residual),1e-10)
        end

        function enstrophyZeroInverseDepthAndFullyNeutralFamilyAreExact(testCase)
            k = 1e-4;
            profiles = {@(z) (5e-3)^2*ones(size(z)),@IMNeutralGeostrophicFiniteDifferenceTests.neutralCap,@(z) zeros(size(z))};
            for i = 1:numel(profiles)
                evp = IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes(N2=profiles{i},zDomain=[-1000 0],k=k,alpha0=9.81*k^2,alphaD=Inf);
                solver = IMSolverFiniteDifference(z=linspace(-1000,0,65).');
                if i==3
                    testCase.verifyError(@() solver.solveEVP(evp,nModes=4),"IMSolver:InsufficientModes")
                    count = 1;
                else
                    count = 4;
                end
                basis = solver.solveEVP(evp,nModes=count);
                zero = find(basis.inverseDepth==0);
                testCase.verifyNumElements(zero,1)
                testCase.verifyEqual(basis.eigenvalues(zero),k^2,AbsTol=0)
                testCase.verifyTrue(isinf(basis.h(zero)))
                G = basis.G([-1000;-500;0]);
                testCase.verifyEqual(G(:,zero),zeros(3,1),AbsTol=1e-10)
                if i==3, testCase.verifyNumElements(basis.h,1); end
            end
            alternate = IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes(N2=@(z) zeros(size(z)), ...
                zDomain=[-1000 0],k=k,alpha0=(9.81*k)*k,alphaD=Inf);
            alternateBasis = IMSolverFiniteDifference(z=linspace(-1000,0,65).').solveEVP(alternate,nModes=1);
            testCase.verifyEqual(alternateBasis.inverseDepth,0,AbsTol=0)
            testCase.verifyTrue(isinf(alternateBasis.h))
            apv = IMInternalModes.geostrophicAPVModes(N2=@(z) zeros(size(z)),zDomain=[-1000 0],g0=0.02,gd=0.03,surfaceBoundary="freeSurface");
            neutral = IMSolverFiniteDifference(z=linspace(-1000,0,33).').solveEVP(apv,nModes=1);
            expected = apv.g*(1/apv.g+1/0.02+1/0.03)/1000;
            testCase.verifyNumElements(neutral.h,1)
            testCase.verifyEqual(neutral.inverseDepth,expected,RelTol=1e-13)
        end

        function zeroAPVUsesOneMultipleRHSFactorizationPerDistinctWavenumber(testCase)
            counter = IMGeostrophicZeroAPVTestCounter();
            problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=@(z) (5e-3)^2*ones(size(z)), ...
                zDomain=[-1000 0],latitude=45,k=[1e-4 2e-4 1e-4 3e-4 2e-4]);
            solver = IMSolverFiniteDifference(z=linspace(-1000,0,65).');
            callback = @(A,B) IMNeutralGeostrophicFiniteDifferenceTests.recordedSolve(counter,A,B);
            basis = IMGeostrophicFiniteDifferenceTools.solveZeroAPV(problem,solver,callback);
            testCase.verifyEqual(counter.callCount,3)
            testCase.verifyEqual(counter.rightHandSideCounts,[2 2 2])
            testCase.verifyEqual(basis.metadata.factorizations,3)
            testCase.verifyEqual(basis.G([-1000;-500;0],pages=1),basis.G([-1000;-500;0],pages=3),AbsTol=0)
        end

        function zeroAPVConstantPairHasSecondOrderShapesAndExactResponses(testCase)
            N0 = 5e-3;
            problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=@(z) N0^2*ones(size(z)),zDomain=[-1000 0],latitude=45,k=[2*pi/1e4 2*pi/1e4]);
            z = linspace(-1000,0,1001).';
            [exactF,exactG] = IMNeutralGeostrophicFiniteDifferenceTests.constantZeroAPV(problem,z,N0);
            errors = zeros(3,2);
            for i = 1:3
                n = 2^(i+5)+1;
                basis = IMSolverFiniteDifference(z=linspace(-1000,0,n).').solveGeostrophicZeroAPVModes(problem);
                F = basis.F(z,pages=1);
                G = basis.G(z,pages=1);
                errors(i,:) = [norm(F-exactF,"fro")/norm(exactF,"fro"),norm(G-exactG,"fro")/norm(exactG,"fro")];
                testCase.verifyEqual(basis.G(0,pages=1)-basis.F(0,pages=1),[1 0],AbsTol=1e-10)
                testCase.verifyEqual(basis.G(-1000,pages=1),[0 1],AbsTol=1e-10)
                testCase.verifyEqual(basis.G(z,pages=1),basis.G(z,pages=2),AbsTol=0)
                testCase.verifyEqual(basis.metadata.factorizations,1)
            end
            testCase.verifyLessThan(max(errors(end,:)),2e-3)
            testCase.verifyGreaterThan(min(errors(1,:)./errors(2,:)),3)
            testCase.verifyGreaterThan(min(errors(2,:)./errors(3,:)),3)
        end

        function zeroAPVNeutralCapConvergesAndKeepsAffineLayer(testCase)
            problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=@IMNeutralGeostrophicFiniteDifferenceTests.neutralCap,zDomain=[-1000 0],latitude=45,k=2*pi/1e4);
            reference = IMSolverSpectral(nEVP=512).solveGeostrophicZeroAPVModes(problem);
            z = linspace(-1000,0,1001).';
            errors = zeros(3,2);
            for i = 1:3
                n = 2^(i+5)+1;
                basis = IMSolverFiniteDifference(z=linspace(-1000,0,n).').solveGeostrophicZeroAPVModes(problem);
                errors(i,:) = [norm(basis.F(z)-reference.F(z),"fro")/norm(reference.F(z),"fro"), ...
                    norm(basis.G(z)-reference.G(z),"fro")/norm(reference.G(z),"fro")];
                cap = linspace(-180,-20,31).';
                F = basis.F(cap);
                G = basis.G(cap);
                testCase.verifyEqual(F,ones(numel(cap),1)*F(1,:),AbsTol=1e-10)
                affine = G(1,:)+(cap-cap(1)).*((G(end,:)-G(1,:))/(cap(end)-cap(1)));
                testCase.verifyEqual(G,affine,AbsTol=1e-10)
            end
            testCase.verifyLessThan(max(errors(end,:)),3e-3)
            testCase.verifyGreaterThan(min(errors(1,:)./errors(2,:)),2.5)
            testCase.verifyGreaterThan(min(errors(2,:)./errors(3,:)),2.5)
        end
    end

    methods (Static, Access = private)
        function values = recordedSolve(counter,A,B)
            counter.record(size(B,2));
            factorization = decomposition(A);
            values = factorization\B;
        end

        function N2 = neutralCap(z)
            t = min(1,max(0,(-z-200)/200));
            N2 = (5e-3)^2*(6*t.^5-15*t.^4+10*t.^3);
        end

        function [F,G] = constantZeroAPV(problem,z,N0)
            % Match independent constant-coefficient exponential solutions.
            k = problem.k(1);
            m = k*N0/abs(problem.f0);
            D = diff(problem.zDomain);
            y = z-problem.zDomain(1);
            pressure = [exp(m*(y-D)),exp(-m*y)];
            displacement = -(problem.g*m/N0^2)*[pressure(:,1),-pressure(:,2)];
            Fs = [1,exp(-m*D)];
            Gs = -(problem.g*m/N0^2)*[1,-exp(-m*D)];
            Gb = -(problem.g*m/N0^2)*[exp(-m*D),-1];
            sigma = double(problem.surfaceBoundary=="freeSurface");
            coefficients = [Gs-sigma*Fs;Gb]\eye(2);
            F = pressure*coefficients;
            G = displacement*coefficients;
        end
    end
end
