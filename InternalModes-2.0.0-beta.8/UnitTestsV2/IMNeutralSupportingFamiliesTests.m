classdef IMNeutralSupportingFamiliesTests < matlab.unittest.TestCase
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
        function fullyNeutralZeroAPVModesHaveExactAffineDisplacement(testCase)
            N2 = @(z) zeros(size(z));
            z = linspace(-1000,0,73).';
            for surfaceBoundary = ["freeSurface","rigidLid"]
                problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=[-1000 0],latitude=45,k=[2*pi/1e4 2*pi/1e3],surfaceBoundary=surfaceBoundary);
                basis = IMSolverSpectral(nEVP=40).solveGeostrophicZeroAPVModes(problem);
                F = basis.F(z);
                G = basis.G(z);
                sigma = double(surfaceBoundary == "freeSurface");
                for iK = 1:numel(problem.k)
                    a = problem.g*problem.k(iK)^2/problem.f0^2;
                    expectedF = [-1 1]/(1000*a+sigma);
                    expectedG = [0 1]-a*(z+1000)*expectedF;
                    testCase.verifyEqual(F(:,:,iK),ones(numel(z),1)*expectedF,AbsTol=1e-12)
                    testCase.verifyEqual(G(:,:,iK),expectedG,AbsTol=1e-10)
                end
                testCase.verifyGreaterThanOrEqual(min(eig(basis.energyMatrix(:,:,1))),-1e-8)
            end
        end

        function neutralCapZeroAPVPairSatisfiesEquationsOffGrid(testCase)
            problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=@IMNeutralSupportingFamiliesTests.neutralCap,zDomain=[-1000 0],latitude=45,k=[2*pi/1e4 2*pi/1e3 2*pi/1e4]);
            counter = IMGeostrophicZeroAPVTestCounter();
            basis = IMInstrumentedGeostrophicZeroAPVSolver(nEVP=1024,solveCounter=counter).solveGeostrophicZeroAPVModes(problem);
            testCase.verifyEqual(counter.callCount,2)
            testCase.verifyEqual(counter.rightHandSideCounts,[2 2])
            z = linspace(-999,-1,501).';
            F = basis.F(z);
            G = basis.G(z);
            h = 0.02;
            Fz = (basis.F(z-2*h)-8*basis.F(z-h)+8*basis.F(z+h)-basis.F(z+2*h))/(12*h);
            Gz = (basis.G(z-2*h)-8*basis.G(z-h)+8*basis.G(z+h)-basis.G(z+2*h))/(12*h);
            for iK = 1:numel(problem.k)
                first = Fz(:,:,iK)+(problem.N2(z)/problem.g).*G(:,:,iK);
                second = Gz(:,:,iK)+(problem.g*problem.k(iK)^2/problem.f0^2)*F(:,:,iK);
                % Physical L2 residuals use off-grid finite differences,
                % independent of the equations used by the solver.
                l2 = @(values) sqrt(trapz(z,sum(abs(values).^2,2)));
                firstScale = l2(Fz(:,:,iK))+l2((problem.N2(z)/problem.g).*G(:,:,iK));
                secondScale = l2(Gz(:,:,iK))+l2((problem.g*problem.k(iK)^2/problem.f0^2)*F(:,:,iK));
                testCase.verifyLessThan(l2(first)/firstScale,1e-7)
                testCase.verifyLessThan(l2(second)/secondScale,1e-7)
                testCase.verifyEqual(basis.G(0,pages=iK)-basis.F(0,pages=iK),[1 0],AbsTol=1e-10)
                testCase.verifyEqual(basis.G(-1000,pages=iK),[0 1],AbsTol=1e-10)
            end
            cap = linspace(-190,-10,17).';
            FCap = basis.F(cap);
            GCap = basis.G(cap);
            for iK = 1:numel(problem.k)
                testCase.verifyLessThan(max(abs(FCap(:,:,iK)-FCap(1,:,iK)),[],"all")/max(abs(F(:,:,iK)),[],"all"),1e-7)
                affine = GCap(1,:,iK)+(cap-cap(1)).*((GCap(end,:,iK)-GCap(1,:,iK))/(cap(end)-cap(1)));
                testCase.verifyEqual(GCap(:,:,iK),affine,AbsTol=1e-7)
            end
            testCase.verifyEqual(basis.F(z,pages=[3 1]),F(:,:,[3 1]),AbsTol=0)
            testCase.verifyEqual(basis.G(z,pages=[3 1]),G(:,:,[3 1]),AbsTol=0)
            rotated = basis.orthogonalizeGeneralizedEnergy(g0=-0.0175,gd=0.0175);
            for iK = 1:numel(problem.k)
                testCase.verifyEqual(rotated.F(z,pages=iK),F(:,:,iK)*rotated.boundaryResponseMatrix(:,:,iK),AbsTol=1e-11)
                testCase.verifyEqual(rotated.G(z,pages=iK),G(:,:,iK)*rotated.boundaryResponseMatrix(:,:,iK),AbsTol=1e-11)
            end
        end

        function fOnlyConstructorRecoversNeutralGWithoutDivision(testCase)
            problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=@(z) zeros(size(z)),zDomain=[-1000 0],latitude=45,k=2*pi/1e4);
            solver = IMSolverSpectral(nEVP=32).configuredForGeostrophicZeroAPVModes(problem);
            a = problem.g*problem.k^2/problem.f0^2;
            values = ones(numel(solver.zNative),1)*([-1 1]/(1000*a+1));
            nativeF = solver.physicalDerivativeMatrix(0)\values;
            compatibility = IMGeostrophicZeroAPVModesBasis(problem=problem,solver=solver,nativeModes=nativeF);
            solved = solver.solveGeostrophicZeroAPVModes(problem);
            z = linspace(-1000,0,31).';
            testCase.verifyEqual(compatibility.G(z),solved.G(z),AbsTol=1e-10)
            testCase.verifyEqual(compatibility.energyMatrix,solved.energyMatrix,RelTol=1e-10,AbsTol=1e-8)
        end

        function neutralCapMDAKeepsConstantModeAndRegularPressure(testCase)
            evp = IMInternalModes.meanDensityAnomalyModes(N2=@IMNeutralSupportingFamiliesTests.neutralCap,zDomain=[-1000 0],g0=0,gd=0);
            basis = IMSolverSpectral(nEVP=384).solveEVP(evp,nModes=6);
            testCase.verifyTrue(any(basis.modeNumber == 0))
            z = linspace(-190,-10,31).';
            testCase.verifyEqual(basis.F(z),zeros(numel(z),6),AbsTol=1e-7)
            G = basis.G(z);
            testCase.verifyLessThan(max(abs(G-G(1,:)),[],"all")/max(abs(G),[],"all"),1e-7)
            testCase.verifyEqual(basis.gramMatrix(variable="G"),diag(basis.signatures),AbsTol=2e-7)
        end

        function neutralSupportedStateRoundTripKeepsComplexPages(testCase)
            N2 = @IMNeutralSupportingFamiliesTests.neutralCap;
            g0 = -0.0175;
            gd = 0.0175;
            evp = IMInternalModes.geostrophicAPVModes(N2=N2,zDomain=[-1000 0],g0=g0,gd=gd,latitude=45,surfaceBoundary="freeSurface");
            basis = IMSolverSpectral(nEVP=384).solveEVP(evp,nModes=8);
            [z,weights] = basis.solver.nativeQuadratureRule([-1000 0]);
            apv = basis.discreteTransform(z=z,weights=weights,nModes=6,variables=["F","G"],gramTolerance=1e-5);
            problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=[-1000 0],latitude=45,k=[2*pi/1e4 2*pi/1e3]);
            zero = IMSolverSpectral(nEVP=384).solveGeostrophicZeroAPVModes(problem);
            transform = IMGeostrophicTransform(apvTransform=apv,zeroAPVModes=zero,g0=g0,gd=gd);
            coefficients = reshape(1:24,6,2,2)*(1+0.3i);
            endpointCoefficients = reshape(1:8,2,2,2)*(0.2-0.7i);
            [state,endpoints] = transform.transformStateBack(APVCoefficients=coefficients,zeroAPVCoefficients=endpointCoefficients);
            [recovered,recoveredEndpoints] = transform.transformStateForward(APV=state,endpointAnomalies=endpoints);
            testCase.verifyEqual(recovered,coefficients,RelTol=1e-9,AbsTol=1e-9)
            testCase.verifyEqual(recoveredEndpoints,endpointCoefficients,RelTol=1e-9,AbsTol=1e-9)
        end

        function zeroAPVRejectsInvalidStratification(testCase)
            profiles = {@(z) -ones(size(z)),@(z) NaN(size(z)),@(z) (1+1i)*ones(size(z))};
            for i = 1:numel(profiles)
                problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=profiles{i},zDomain=[-1 0],k=1);
                testCase.verifyError(@() IMSolverSpectral(nEVP=32).solveGeostrophicZeroAPVModes(problem),"IMGeostrophicZeroAPVModes:InvalidStratification")
            end
        end
    end

    methods (Static, Access = private)
        function N2 = neutralCap(z)
            t = min(1,max(0,(-z-200)/200));
            N2 = (5e-3)^2*(6*t.^5-15*t.^4+10*t.^3);
        end
    end
end
