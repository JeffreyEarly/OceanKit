classdef IMBasisCollectionTests < matlab.unittest.TestCase
    properties (Access = private)
        originalPath
    end
    methods (TestClassSetup)
        function addRepositoryPath(testCase)
            testCase.originalPath = path;
            addpath(fileparts(fileparts(mfilename("fullpath"))));
        end
    end
    methods (TestClassTeardown)
        function restoreRepositoryPath(testCase)
            path(testCase.originalPath);
        end
    end
    methods (Test)
        function scalarNormalizationAndScientificLabelsAreSnapshots(testCase)
            basis = testCase.scalarBasis(3);
            z = linspace(-1,0,15).';
            collection = IMBasisCollection({basis},basisIndex=[1 1]);
            basis = basis.addNormalization("twice",@(b,j) 2*b.innerProductNormFactor(j));
            basis.normalization = "twice";
            expected = 2*basis.u(z);
            values = collection.evaluate(z,variable="u");
            testCase.verifyEqual(values(:,:,1),expected,AbsTol=1e-13);
            testCase.verifyEqual(values(:,:,2),expected,AbsTol=1e-13);
            derivative = collection.evaluate(z,variable="u",derivativeOrder=1,pages=2,columns=[3 1]);
            expectedDerivative = 2*basis.uz(z);
            testCase.verifyEqual(derivative,expectedDerivative(:,[3 1]),AbsTol=1e-13);
            testCase.verifyEqual(collection.metadata.normalization,"rms");
            testCase.verifyEqual(collection.metadata.columnLabels,basis.columnLabels);
            testCase.verifyEmpty(collection.kappa);
            testCase.verifyError(@() collection.evaluate(z,variable="u",derivativeOrder=2),"IMBasisCollection:UnsupportedDerivative");
            testCase.verifyError(@() collection.evaluate(z,variable="F"),"IMBasisCollection:UnsupportedVariable");
        end

        function analyticalAndNumericalModesShareThePageConvention(testCase)
            N0 = 5e-3;
            zDomain = [-1000 0];
            N2 = @(z) N0^2*ones(size(z));
            evp = IMInternalModes.hydrostaticGModes(rotationRate=0,N2=N2,zDomain=zDomain,surfaceBoundary=IMBoundaryCondition.dirichlet());
            numerical = IMSolverSpectral(nEVP=48).solveEVP(evp,nModes=3);
            solution = IMConstantStratificationSolution(rotationRate=0,N0=N0,zDomain=zDomain);
            analytical = solution.internalModes(evp,nModes=3);
            collection = IMBasisCollection({numerical,analytical},basisIndex=[2 1 2]);
            z = linspace(-1000,0,19).';
            values = collection.evaluate(z,variable="G",pages=[2 3 1]);
            testCase.verifyEqual(values(:,:,1),values(:,:,2),RelTol=1e-10,AbsTol=2e-10);
            testCase.verifyEqual(values(:,:,2),values(:,:,3),AbsTol=0);
            derivatives = collection.evaluate(z,variable="G",derivativeOrder=1);
            testCase.verifyEqual(derivatives(:,:,1),derivatives(:,:,2),AbsTol=1e-11);
            testCase.verifyEqual(collection.metadata(1).variables,["F" "G"]);
            testCase.verifyEqual(collection.metadata(2).derivativeOrders,{[0 1],0:2});
            derivativeF = collection.evaluate(z,variable="F",derivativeOrder=1);
            testCase.verifyEqual(derivativeF(:,:,1),numerical.evaluate(z,variable="F",derivativeOrder=1),AbsTol=1e-13);
            testCase.verifyEqual(derivativeF(:,:,2),analytical.evaluate(z,variable="F",derivativeOrder=1),AbsTol=1e-13);
        end

        function waveKappaMappingsRejectApproximateSubstitution(testCase)
            N2 = @(z) 2.5e-5*ones(size(z));
            solver = IMSolverSpectral(nEVP=40);
            basis1 = solver.solveEVP(IMInternalModes.waveModesAtWavenumber(rotationRate=0,N2=N2,zDomain=[-1000 0],k=1e-4,surfaceBoundary=IMBoundaryCondition.dirichlet()),nModes=3);
            basis2 = solver.solveEVP(IMInternalModes.waveModesAtWavenumber(rotationRate=0,N2=N2,zDomain=[-1000 0],k=2e-4,surfaceBoundary=IMBoundaryCondition.dirichlet()),nModes=3);
            collection = IMBasisCollection({basis1,basis2},kappa=[2e-4 1e-4 2e-4],basisIndex=[2 1 2]);
            testCase.verifyEqual(collection.kappa,[2e-4 1e-4 2e-4]);
            values = collection.evaluate([-700;-100],variable="F");
            testCase.verifyEqual(values(:,:,1),values(:,:,3),AbsTol=0);
            testCase.verifyEqual(values(:,:,2),basis1.F([-700;-100]),AbsTol=0);
            inferred = IMBasisCollection({basis2,basis1});
            testCase.verifyEqual(inferred.kappa,[2e-4 1e-4]);
            testCase.verifyError(@() IMBasisCollection({basis1},kappa=1.0001e-4),"IMBasisCollection:KappaMismatch");
        end

        function boundaryPagesPreserveEndpointResponsesAndRotations(testCase)
            N0 = 5e-3;
            solution = IMConstantStratificationSolution(N0=N0,zDomain=[-1000 0],rotationRate=1e-4/(2*sind(33)),latitude=33);
            analytical = solution.geostrophicZeroAPVModesAtWavenumber([1e-4 2e-4]);
            problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=@(z) N0^2*ones(size(z)),zDomain=[-1000 0],rotationRate=1e-4/(2*sind(33)),latitude=33,k=[1e-4 2e-4]);
            numerical = IMSolverSpectral(nEVP=48).solveGeostrophicZeroAPVModes(problem);
            for basis = {analytical,numerical}
                collection = IMBasisCollection(basis,basisIndex=[1 1 1],sourcePage=[2 1 2]);
                testCase.verifyEqual(collection.kappa,[2e-4 1e-4 2e-4]);
                testCase.verifyEqual(collection.metadata.columnKind,"endpoint");
                testCase.verifyEqual(collection.metadata.endpointLabels,["surface" "bottom"]);
                testCase.verifyEmpty(collection.metadata.modeNumber);
                F = collection.evaluate([0;-1000],variable="F");
                G = collection.evaluate([0;-1000],variable="G");
                response = G;
                response(1,:,:) = G(1,:,:)-F(1,:,:);
                testCase.verifyEqual(response,repmat(eye(2),1,1,3),AbsTol=1e-9);
                derivative = collection.evaluate([0;-1000],variable="F",derivativeOrder=1);
                testCase.verifyEqual(derivative,-N0^2/basis{1}.g*G,AbsTol=1e-15);
                selected = collection.evaluate([0;-1000],variable="G",pages=[3 1],columns=2);
                testCase.verifyEqual(selected,G(:,2,[3 1]),AbsTol=0);
                expectedDerivative = basis{1}.evaluate(0,variable="G",derivativeOrder=1,pages=[2 1 2]);
                testCase.verifyEqual(collection.evaluate(0,variable="G",derivativeOrder=1),expectedDerivative,AbsTol=1e-12);
            end
            rotated = analytical.orthogonalizeGeneralizedEnergy(g0=0,gd=0);
            rotatedCollection = IMBasisCollection({rotated});
            testCase.verifyEqual(rotatedCollection.metadata.columnLabels,rotated.columnLabels);
            testCase.verifyEqual(rotatedCollection.metadata.endpointLabels,["surface" "bottom"]);
            testCase.verifyEqual(rotatedCollection.metadata.boundaryResponseMatrixStatus,"availableOnBasis");
            testCase.verifyEqual(rotatedCollection.bases{1}.boundaryResponseMatrix,rotated.boundaryResponseMatrix);
        end

        function meanAndAPVBasesRetainTheirFamilyIdentity(testCase)
            N2 = @(z) 2.5e-5*ones(size(z));
            solver = IMSolverSpectral(nEVP=48);
            apv = solver.solveEVP(IMInternalModes.geostrophicAPVModes(rotationRate=0,N2=N2,zDomain=[-1000 0],g0=0,gd=0,surfaceBoundary="freeSurface"),nModes=3);
            mean = solver.solveEVP(IMInternalModes.meanDensityAnomalyModes(rotationRate=0,N2=N2,zDomain=[-1000 0],g0=0,gd=0),nModes=3);
            for basis = {apv,mean}
                collection = IMBasisCollection(basis);
                variable = string(basis{1}.evp.formulation);
                testCase.verifyEqual(collection.metadata.family,string(basis{1}.evp.name));
                testCase.verifyEqual(collection.metadata.modeNumber,basis{1}.modeNumber);
                testCase.verifyEqual(collection.evaluate([-700;-100],variable=variable,derivativeOrder=1),basis{1}.uz([-700;-100]),AbsTol=0);
            end
        end

        function assessmentsKeepUnequalPerPageCounts(testCase)
            collection = IMBasisCollection({testCase.scalarBasis(4),testCase.scalarBasis(2)},basisIndex=[2 1 2]);
            z = linspace(-1,0,17).';
            weights = [0.5;ones(15,1);0.5]/16;
            smaller = collection.assess(z,weights,page=1,prefixColumnCounts=1:2);
            larger = collection.assess(z,weights,page=2,prefixColumnCounts=[2 4]);
            testCase.verifyClass(smaller,"struct");
            testCase.verifyEqual(smaller.projection.columnCount,2);
            testCase.verifyEqual(larger.projection.columnCount,4);
            testCase.verifyEqual(smaller.identity.basisIndex,2);
            testCase.verifyEqual(larger.identity.basisIndex,1);
            checkBasisAssessment(smaller,gramTolerance=0);
            testCase.verifyEqual(collection.projectionOnGrid(z,weights,variable="u",page=2).columnCount,4);
            testCase.verifyEqual(collection.projectionOnGrid(z,weights,variable="u",page=3,prepareRecovery=true).recovery.forwardMatrix,smaller.projection.recovery.forwardMatrix);
        end

        function invalidSelectionsAndHeterogeneousCountsAreExplicit(testCase)
            basis3 = testCase.scalarBasis(3);
            basis2 = testCase.scalarBasis(2);
            collection = IMBasisCollection({basis3,basis2});
            testCase.verifyError(@() collection.evaluate(-0.5,variable="u"),"IMBasisCollection:HeterogeneousColumns");
            testCase.verifySize(collection.evaluate(-0.5,variable="u",columns=[2 1]),[1 2 2]);
            testCase.verifyError(@() collection.evaluate(-0.5,variable="u",columns=3),"IMBasisCollection:InvalidColumn");
            testCase.verifyError(@() collection.evaluate(-0.5,variable="u",pages=3),"IMBasisCollection:InvalidPage");
            testCase.verifyError(@() collection.evaluate(1,variable="u",pages=1),"IMBasisCollection:OutsideDomain");
            testCase.verifyError(@() IMBasisCollection({}),"IMBasisCollection:EmptyCollection");
            testCase.verifyError(@() IMBasisCollection({1}),"IMBasisCollection:UnsupportedBasis");
            testCase.verifyError(@() IMBasisCollection({basis3},basisIndex=2),"IMBasisCollection:InvalidBasisIndex");
            testCase.verifyError(@() IMBasisCollection({basis3},sourcePage=2),"IMBasisCollection:InvalidSourcePage");
            testCase.verifyError(@() IMBasisCollection({basis3},kappa=[1 2]),"IMBasisCollection:InvalidKappaCount");
            duplicateLabels = IMBasisSet(solver=basis3.solver,evp=basis3.evp,nativeModes=basis3.nativeModes,eigenvalues=basis3.eigenvalues,modeNumber=[1 1 3]);
            testCase.verifyError(@() IMBasisCollection({duplicateLabels}),"IMBasisCollection:InvalidColumnLabels");
        end
    end
    methods (Static, Access = private)
        function basis = scalarBasis(nModes)
            evp = IMEigenvalueProblem(zDomain=[-1 0],p=1,q=0,r=1,surfaceBoundary=IMBoundaryCondition.dirichlet(),bottomBoundary=IMBoundaryCondition.dirichlet());
            basis = IMSolverSpectral(nEVP=32).solveEVP(evp,nModes=nModes);
        end
    end
end
