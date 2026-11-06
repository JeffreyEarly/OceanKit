classdef IMHydrostaticNativeBasisTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addProvider(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
        end
    end

    methods (Test)
        function analyticalSectorRetainsNativeMetricsAndDerivatives(testCase)
            N0 = 5e-3; D = 1000;
            source = testCase.source(N0,D);
            basis = IMHydrostaticModesBasis(source,normalization="rmsF");
            z = linspace(-D,0,33).';
            F = basis.F(z); G = basis.G(z);
            testCase.verifyClass(basis,"IMHydrostaticModesBasis");
            testCase.verifyEqual(basis.nColumns,5);
            testCase.verifyEqual(basis.modeNumber,0:4);
            testCase.verifyEqual(basis.h,[Inf source.h]);
            testCase.verifyEqual(F(:,1),ones(size(z)),AbsTol=1e-14);
            testCase.verifyEqual(G(:,1),zeros(size(z)));
            testCase.verifyEqual(G([1 end],:),zeros(2,5));
            derivative = basis.Gz(z);
            testCase.verifyEqual(F(:,2:end),basis.h(2:end).*derivative(:,2:end),AbsTol=1e-13);
            testCase.verifyEqual(basis.evaluate(z,variable="F",derivativeOrder=1),-N0^2/source.evp.g*G,AbsTol=1e-12);
            testCase.verifyEqual(basis.gramMatrix(variable="F"),D*eye(5),AbsTol=1e-8);
            testCase.verifyEqual(basis.gramMatrix(variable="G"),diag([0 D./source.h]),AbsTol=1e-5);
            selected = basis.selectColumns([4 1 2]);
            testCase.verifyClass(selected,"IMHydrostaticModesBasis");
            testCase.verifyEqual(selected.modeNumber,[3 0 1]);
            testCase.verifyEqual(selected.F(z),F(:,[4 1 2]),AbsTol=1e-13);
            testCase.verifyEqual(selected.G(z),G(:,[4 1 2]),AbsTol=1e-10);
        end

        function fixedGridPreservesBasisIdentityAndBothActiveChannels(testCase)
            N0 = 5e-3; D = 1000;
            basis = IMHydrostaticModesBasis(testCase.source(N0,D));
            grid = IMVerticalGrid.physicalDepth([-D 0],65);
            transform = basis.discreteTransform(grid);
            testCase.verifyClass(transform.basis,"IMHydrostaticModesBasis");
            testCase.verifyEqual(transform.modeNumber,basis.modeNumber);
            testCase.verifyEqual(transform.activeModeMask(variable="G"),[false true(1,4)]);
            testCase.verifyEqual(transform.inverseMatrix(variable="F"),basis.F(grid.z),AbsTol=0);
            testCase.verifyEqual(transform.inverseMatrix(variable="G"),basis.G(grid.z),AbsTol=0);
            coefficient = [2;-1;0.5;0.2;-0.1];
            testCase.verifyEqual(transform.transformForward(basis.F(grid.z)*coefficient,variable="F"),coefficient,AbsTol=1e-12);
            expectedG = coefficient; expectedG(1) = 0;
            testCase.verifyEqual(transform.transformForward(basis.G(grid.z)*coefficient,variable="G"),expectedG,AbsTol=1e-12);
            collection = IMBasisCollection({basis},basisIndex=[1 1]);
            testCase.verifyEqual(collection.metadata.derivativeOrders,{[0 1],0:2});
            expected = basis.evaluate(grid.z,variable="F",derivativeOrder=1);
            actual = collection.evaluate(grid.z,variable="F",derivativeOrder=1,sampleChunkSize=8);
            testCase.verifyEqual(actual,repmat(expected,1,1,2),AbsTol=1e-15);
        end

        function auxiliaryLeakageSolveRetainsTheBarotropicColumn(testCase)
            N0 = 5e-3; D = 1000;
            basis = IMHydrostaticModesBasis(testCase.source(N0,D));
            internal = basis.solver.solveEVP(basis.evp,nModes=8);
            extended = IMHydrostaticModesBasis(internal);
            z = linspace(-D,0,65).';
            weights = D/64*[0.5;ones(63,1);0.5];
            [transform,assessment] = basis.discreteTransform(z=z,weights=weights,nModes=3,nCheckModes=8,leakageTolerance=1e-8);
            [reference,referenceAssessment] = extended.discreteTransform(z=z,weights=weights,nModes=3,nCheckModes=8,leakageTolerance=1e-8);
            testCase.verifyEqual(transform.modeNumber,0:2);
            testCase.verifyEqual(transform.inverseMatrix(variable="F"),reference.inverseMatrix(variable="F"),AbsTol=1e-12);
            testCase.verifyEqual(assessment.leakagePolicy.error,referenceAssessment.leakagePolicy.error,AbsTol=1e-12);
            testCase.verifyTrue(all(assessment.leakagePolicy.accepted));
        end
    end

    methods (Static, Access = private)
        function basis = source(N0,D)
            evp = IMInternalModes.hydrostaticGModes(N2=@(z) N0^2*ones(size(z)),zDomain=[-D 0],rotationRate=0,surfaceBoundary=IMBoundaryCondition.dirichlet());
            basis = IMSolverSpectral(nEVP=48,coordinateKind="z").solveEVP(evp,nModes=4);
        end
    end
end
