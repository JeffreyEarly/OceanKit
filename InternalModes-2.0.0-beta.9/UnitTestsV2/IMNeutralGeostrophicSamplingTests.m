classdef IMNeutralGeostrophicSamplingTests < matlab.unittest.TestCase
    properties (Access = private)
        originalPath
    end

    methods (TestClassSetup)
        function addRepositoryPath(testCase)
            repoRoot = fileparts(fileparts(mfilename("fullpath")));
            testCase.originalPath = path;
            addpath(repoRoot,fullfile(repoRoot,"UnitTestsV2","TestSupport"));
        end
    end

    methods (TestClassTeardown)
        function restoreRepositoryPath(testCase)
            path(testCase.originalPath);
        end
    end

    methods (Test)
        function mdaSamplingRetainsTheHighestPressureCoefficient(testCase)
            for n = [12 17]
                evp = IMInternalModes.meanDensityAnomalyModes(N2=@(z) 9.81*ones(size(z)), ...
                    zDomain=[-2 0],g0=Inf,gd=Inf);
                solver = IMSolverSpectral(nEVP=n,coordinateKind="z").configuredForEVP(evp);
                G = zeros(n,1);
                G(end) = 1;
                basis = IMMeanDensityAnomalyModesBasis(solver=solver,evp=evp,nativeModes=G,eigenvalues=3,modeNumber=1);
                basis = basis.addNormalization("raw",@(~,~) 1);
                basis.normalization = "raw";
                % Integral of T_(n-1): T_n/(2*n)-T_(n-2)/(2*(n-2)).
                primitive = zeros(n+1,1);
                primitive(end) = 1/(2*n);
                primitive(end-2) = -1/(2*(n-2));
                F = -primitive;
                F(1) = sum(primitive);
                z = linspace(-2,0,101).';
                sampled = sampleNeutralGeostrophicBasis(basis,z);
                expectedF = solver.evaluateNativeModes(F,z);
                expectedFz = solver.evaluatePhysicalDerivative(F,z,1);
                expectedGz = solver.evaluatePhysicalDerivative(G,z,1);
                expectedGzz = solver.evaluatePhysicalDerivative(G,z,2);

                testCase.verifyEqual(sampled.F,expectedF,AbsTol=2e-14)
                testCase.verifyEqual(sampled.Fz,expectedFz,AbsTol=2e-12)
                testCase.verifyEqual(sampled.Gz,expectedGz,AbsTol=2e-12)
                testCase.verifyEqual(sampled.Gzz,expectedGzz,AbsTol=2e-11)
                testCase.verifyGreaterThan(max(abs(sampled.Gz-3*sampled.F)),1)
            end
        end

        function pairedSamplingDifferentiatesStoredFieldsIndependently(testCase)
            n = 16;
            evp = IMInternalModes.geostrophicAPVModes(N2=@(z) ones(size(z)),zDomain=[-2 0],g0=Inf,gd=Inf,surfaceBoundary="freeSurface");
            solver = IMSolverSpectral(nEVP=n,coordinateKind="z").configuredForEVP(evp);
            F = zeros(n,1);
            F([1 4]) = [1 .2];
            G = zeros(n,1);
            G([3 5]) = [-.5 2];
            basis = IMGeostrophicModesBasis(solver=solver,evp=evp,nativeModes=F,nativeGModes=G, ...
                eigenvalues=7,inverseDepth=7,modeNumber=1);
            basis = basis.addNormalization("raw",@(~,~) 1);
            basis.normalization = "raw";
            z = linspace(-2,0,73).';
            sampled = sampleNeutralGeostrophicBasis(basis,z);

            testCase.verifyEqual(sampled.Fz,solver.evaluatePhysicalDerivative(F,z,1),AbsTol=2e-13)
            testCase.verifyEqual(sampled.Gz,solver.evaluatePhysicalDerivative(G,z,1),AbsTol=2e-13)
            testCase.verifyGreaterThan(max(abs(sampled.Gz-7*sampled.F)),1)
        end

        function legacyZeroAPVSamplingRetainsTheIntegratedCompanionDegree(testCase)
            n = 12;
            problem = IMGeostrophicZeroAPVModes(N2=@(z) ones(size(z)),zDomain=[-2 0], ...
                k=1,rotationRate=.5,latitude=90,endpoints="surface");
            solver = IMSolverSpectral(nEVP=n,coordinateKind="z").configuredForGeostrophicZeroAPVModes(problem);
            F = zeros(n,1);
            F(end) = 1;
            basis = IMGeostrophicZeroAPVModesBasis(problem=problem,solver=solver,nativeModes=F);
            primitive = zeros(n+1,1);
            primitive(end) = 1/(2*n);
            primitive(end-2) = -1/(2*(n-2));
            G = -problem.g*primitive;
            G(1) = problem.g*sum(primitive)+2;
            z = linspace(-2,0,101).';
            sampled = sampleNeutralGeostrophicBasis(basis,z);

            testCase.verifyEqual(sampled.G,solver.evaluateNativeModes(G,z),AbsTol=2e-13)
            testCase.verifyEqual(sampled.Gz,solver.evaluatePhysicalDerivative(G,z,1),AbsTol=2e-11)
            testCase.verifyEqual(sampled.Fz,solver.evaluatePhysicalDerivative(F,z,1),AbsTol=2e-11)
        end

        function piecewiseStructuresPreserveTheirStoredDerivativeFields(testCase)
            basis = struct("evaluate",@IMNeutralGeostrophicSamplingTests.syntheticFields);
            z = [-2;-1.3;-.2;0];
            expected = IMNeutralGeostrophicSamplingTests.syntheticFields(z);
            sampled = sampleNeutralGeostrophicBasis(basis,z);
            expected.originalShape = size(expected.F);

            testCase.verifyEqual(sampled,expected)
        end
    end

    methods (Static,Access = private)
        function values = syntheticFields(z)
            values = struct("F",[z.^2 3*z+1],"G",[cos(z) exp(z)], ...
                "Fz",[2*z 3*ones(size(z))],"Gz",[-sin(z) exp(z)], ...
                "Gzz",[-cos(z) exp(z)],"extraField","preserved");
        end
    end
end
