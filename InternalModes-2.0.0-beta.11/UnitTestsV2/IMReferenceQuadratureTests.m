classdef IMReferenceQuadratureTests < matlab.unittest.TestCase
    % Compare public profile normalizations with independent continuous integrals.
    methods (TestClassSetup)
        function configurePath(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
        end
    end

    methods (Test)
        function finiteDifferenceScalarRmsUsesEvaluatedProfiles(testCase)
            zNative = linspace(-1,0,24).';
            solver = IMSolverFiniteDifference(z=zNative);
            evp = IMEigenvalueProblem(zDomain=[-1 0],p=1,q=0,r=1);
            basis = solver.solveEVP(evp,nModes=22);
            [z,weights] = fourPointCellRule(zNative);
            measuredRms = sqrt(sum(weights.*abs(basis.u(z)).^2,1));
            testCase.verifyEqual(measuredRms,ones(1,22),AbsTol=2e-12);
            selected = basis.selectColumns([22 20]);
            values = basis.u(z);
            testCase.verifyEqual(selected.u(z),values(:,[22 20]),AbsTol=2e-12);
            % Native quadrature continues to describe the nodal discretization.
            [nativePoints,nativeWeights] = solver.nativeQuadratureRule([-1 0]);
            testCase.verifyEqual(nativePoints,zNative);
            testCase.verifyEqual(nativeWeights,[.5;ones(22,1);.5]/23,AbsTol=1e-15);
        end

        function finiteDifferenceWaveNormalizationsUseBothEvaluatedProfiles(testCase)
            zNative = linspace(-1000,0,24).';
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            k = 1e-3;
            basis = IMSolverFiniteDifference(z=zNative).solveEVP(IMInternalModes.waveModesAtWavenumber(background,k=k),nModes=22);
            [z,weights] = fourPointCellRule(zNative);
            for name = ["rmsF","rmsG","waveEnergy"]
                normalized = basis.withNormalization(name);
                F = normalized.F(z);
                G = normalized.G(z);
                switch name
                    case "rmsF", integrand = abs(F).^2;
                    case "rmsG", integrand = abs(G).^2;
                    case "waveEnergy", integrand = abs(F).^2+(k*1000)^2*abs(G).^2;
                end
                testCase.verifyEqual(sum(weights.*integrand,1)/1000,ones(1,22),AbsTol=3e-12);
            end
        end

        function staggeredGeostrophicRmsUsesBothProfileKnotSets(testCase)
            zNative = linspace(-1000,0,24).';
            evp = IMInternalModes.geostrophicAPVModes(N2=@(z) 1e-4*ones(size(z)),zDomain=[-1000 0],g0=0,gd=0,surfaceBoundary="freeSurface");
            basis = IMSolverFiniteDifference(z=zNative).solveEVP(evp,nModes=20);
            knots = sort([zNative;(zNative(1:end-1)+zNative(2:end))/2]);
            [z,weights] = fourPointCellRule(knots);
            for variable = ["F","G"]
                normalized = basis.withNormalization("rms"+variable);
                values = normalized.evaluate(z,variable=variable);
                testCase.verifyEqual(sum(weights.*abs(values).^2,1)/1000,ones(1,20),AbsTol=3e-12);
            end
        end

        function spectralPartialGramUsesTheOriginalCoordinateJacobian(testCase)
            background = IMBackground.exponentialStratification(N0=.01,b=700,zDomain=[-1000 0],surfaceCondition="rigidLid");
            basis = IMSolverSpectral(nEVP=48,coordinateKind="wkb").solveEVP(IMInternalModes.hydrostaticGModes(background),nModes=4);
            bounds = [-820 -130];
            expected = integral(@(z) basis.F(z).'*basis.F(z),bounds(1),bounds(2),ArrayValued=true,AbsTol=1e-10,RelTol=1e-12);
            testCase.verifyEqual(basis.gramMatrix(variable="F",zBounds=bounds),expected,AbsTol=3e-9);
        end
    end
end

function [z,weights] = fourPointCellRule(knots)
% Independent closed-form Gauss rule, exact for products of cubic profiles.
outer = sqrt((3+2*sqrt(6/5))/7);
inner = sqrt((3-2*sqrt(6/5))/7);
nodes = [-outer;-inner;inner;outer];
rule = [(18-sqrt(30))/36;(18+sqrt(30))/36;(18+sqrt(30))/36;(18-sqrt(30))/36];
widths = diff(knots).';
z = knots(1:end-1).'+((nodes+1)/2).*widths;
weights = (rule/2).*widths;
z = z(:);
weights = weights(:);
end
