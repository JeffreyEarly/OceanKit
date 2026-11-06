classdef IMProjectionOnGridTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function configurePath(testCase)
            originalPath = path;
            addpath(fileparts(fileparts(mfilename("fullpath"))));
            testCase.addTeardown(@() path(originalPath));
        end
    end
    methods (Test)
        function physicalProjectionIsTheDefaultAndRecoveryIsExplicit(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            basis = IMSolverAnalytical().solveEVP(IMInternalModes.hydrostaticFModes(background),nModes=4);
            z = linspace(-1000,0,19).';
            weights = 1000/18*[.5;ones(17,1);.5];
            for source = {basis,IMBasisCollection({basis})}
                physical = source{1}.projectionOnGrid(z,weights,variable="F");
                testCase.verifyEmpty(physical.recovery);
                testCase.verifyFalse(isprop(physical,"forwardMatrix"));
                testCase.verifyFalse(ismethod(physical,"project"));
                testCase.verifySize(physical.projectPhysical(exp(z/200)),[4 1]);
                prepared = source{1}.projectionOnGrid(z,weights,variable="F",prepareRecovery=true);
                testCase.verifyTrue(prepared.recovery.success);
                testCase.verifyEqual(prepared.physicalDualMatrix,physical.physicalDualMatrix);
                testCase.verifyEqual(prepared.recovery.transformForward(basis.F(z)),eye(4),AbsTol=1e-12);
            end
        end

        function collectionsKeepAuxiliaryWavePairingsInAssessment(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            problems = {IMInternalModes.waveModesAtWavenumber(background,k=0),IMInternalModes.waveModesAtWavenumber(background,k=1e-3),IMInternalModes.waveModesAtFrequency(background,omega=.005)};
            z = linspace(-1000,0,65).';
            weights = 1000/64*[.5;ones(63,1);.5];
            for solver = {IMSolverSpectral(nEVP=48),IMSolverAnalytical()}
                for problem = problems
                    basis = solver{1}.solveEVP(problem{1},nModes=4);
                    collection = IMBasisCollection({basis});
                    for source = {basis,collection}
                        testCase.verifyError(@() source{1}.projectionOnGrid(z,weights,variable="G"),"IMScalarProjection:UnavailableProjection");
                        testCase.verifyError(@() source{1}.projectionOnGrid(z,weights,variable="G",prepareRecovery=true),"IMScalarProjection:UnavailableProjection");
                    end
                    assessment = collection.assess(z,weights,variable="G");
                    testCase.verifyEqual(assessment.identity.variable,"G");
                    testCase.verifyTrue(assessment.projection.recovery.success);
                    testCase.verifyEqual(assessment.projection.recovery.transformForward(basis.G(z)),eye(4),AbsTol=1e-11);
                end
            end
        end

        function collectionPhysicalSelectionMatchesRequestedBasis(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            first = IMSolverAnalytical().solveEVP(IMInternalModes.hydrostaticGModes(background),nModes=4);
            second = first.selectColumns([4 2 1]);
            collection = IMBasisCollection({first,second},basisIndex=[2 1 2]);
            z = linspace(-1000,0,33).';
            weights = 1000/32*[.5;ones(31,1);.5];
            for page = 1:3
                basis = collection.bases{collection.basisIndex(page)};
                for columns = {":",[3 1],[]}
                    indices = columns{1};
                    if isequal(indices,":"), indices = 1:basis.nColumns; end
                    expected = basis.projectionOnGrid(z,weights,variable="G",columns=indices);
                    actual = collection.projectionOnGrid(z,weights,page=page,variable="G",columns=columns{1});
                    testCase.verifyEqual(actual.physicalDualMatrix,expected.physicalDualMatrix);
                    testCase.verifyEqual(actual.columnLabels,expected.columnLabels);
                    testCase.verifyEmpty(actual.recovery);
                end
            end
            testCase.verifyError(@() collection.projectionOnGrid(z,weights,variable="G",page=4),"IMBasisCollection:InvalidPage");
            testCase.verifyError(@() collection.projectionOnGrid(z,weights,variable="G",columns=[1;2]),"IMProjection:InvalidColumns");
        end

        function scalarAndAlignedRecipesPreserveExistingTransforms(testCase)
            solver = IMSolverSpectral(nEVP=48);
            scalarEVP = IMEigenvalueProblem(zDomain=[-1 0],p=1,q=0,r=1,surfaceBoundary=IMBoundaryCondition.dirichlet(),bottomBoundary=IMBoundaryCondition.dirichlet());
            scalar = solver.solveEVP(scalarEVP,nModes=4);
            z = linspace(-1,0,25).';
            weights = [0.5;ones(23,1);0.5]/24;
            old = scalar.discreteTransform(z=z,weights=weights,nModes=3);
            collection = IMBasisCollection({scalar});
            projection = collection.projectionOnGrid(z,weights,columns=1:3,variable="u",prepareRecovery=true);
            testCase.verifyEqual(projection.sampledBasis,old.inverseMatrix,AbsTol=1e-12)
            testCase.verifyEqual(projection.metricMatrix,old.metricMatrix,AbsTol=1e-12)
            testCase.verifyEqual(projection.targetGramMatrix,old.targetGramMatrix,AbsTol=1e-12)
            testCase.verifyEqual(projection.recovery.forwardMatrix,old.forwardMatrix,AbsTol=1e-11)
            evp = IMInternalModes.hydrostaticFModes(rotationRate=0,N2=@(z) ones(size(z)),zDomain=[-1 0],surfaceBoundary=IMBoundaryCondition.neumann());
            aligned = solver.solveEVP(evp,nModes=4);
            old = aligned.discreteTransform(z=z,weights=weights,nModes=4,gramTolerance=1);
            for variable = ["F","G"]
                projection = IMBasisCollection({aligned}).projectionOnGrid(z,weights,variable=variable,prepareRecovery=true);
                testCase.verifyEqual(projection.sampledBasis,old.inverseMatrix(variable=variable),AbsTol=1e-12)
                testCase.verifyEqual(projection.recovery.forwardMatrix,old.forwardMatrix(variable=variable),AbsTol=1e-10)
                testCase.verifyEqual(projection.activeColumnMask,old.activeModeMask(variable=variable))
            end
        end
        function analyticalRecipesRetainNormalizationAndReferenceProvenance(testCase)
            evp = IMInternalModes.hydrostaticGModes(rotationRate=0,N2=@(z) 1e-4*ones(size(z)),zDomain=[-1000 0],surfaceBoundary=IMBoundaryCondition.dirichlet());
            basis = IMConstantStratificationSolution(rotationRate=0,N0=0.01,zDomain=[-1000 0]).internalModes(evp,nModes=4);
            collection = IMBasisCollection({basis});
            z = linspace(-1000,0,21).';
            weights = [0.5;ones(19,1);0.5]*50;
            projection = collection.projectionOnGrid(z,weights,variable="G",prepareRecovery=true);
            expected = basis.G(z);
            basis.normalization = "maxAbsG";
            testCase.verifyEqual(collection.evaluate(z,variable="G"),expected,AbsTol=1e-12)
            testCase.verifyEqual(collection.evaluate(z,variable="G",columns=[3 1]),expected(:,[3 1]),AbsTol=1e-12)
            testCase.verifyEqual(projection.provenance.referencePointCount,128)
            testCase.verifyEqual(projection.provenance.referenceIntegration,"Gauss-Legendre rule")
            testCase.verifyEqual(projection.provenance.referenceConvergence,"unverified")
            testCase.verifyEqual(projection.provenance.solveAccuracy,"unverified")
            testCase.verifyEqual(projection.prefix(2).provenance,projection.provenance)
            testCase.verifyEqual(collection.assess(z,weights,variable="G").identity.projectionProvenance,projection.provenance)
            testCase.verifyError(@() collection.projectionOnGrid(z,weights,variable="G",columns=[1 1]),"IMAnalyticalInternalModesBasis:InvalidColumns")
            largerBasis = IMConstantStratificationSolution(rotationRate=0,N0=0.01,zDomain=[-1000 0]).internalModes(evp,nModes=40);
            largerRecipe = largerBasis.projectionRecipe(variable="G");
            selectedRecipe = largerBasis.selectColumns([1 3]).projectionRecipe(variable="G");
            testCase.verifyGreaterThan(largerRecipe.provenance.referencePointCount,128)
            testCase.verifyEqual(selectedRecipe.provenance,largerRecipe.provenance)
            testCase.verifyEqual(selectedRecipe.targetGramMatrix,largerRecipe.targetGramMatrix([1 3],[1 3]),AbsTol=1e-10)
        end
        function signedEndpointMetricsAndMissingTracesAreExplicit(testCase)
            solver = IMSolverSpectral(nEVP=48);
            evp = IMEigenvalueProblem(zDomain=[-1 0],p=1,q=0,r=1,surfaceBoundary=IMBoundaryCondition(a=0,b=1,c=1,d=0),bottomBoundary=IMBoundaryCondition.dirichlet());
            basis = solver.solveEVP(evp,nModes=3);
            z = linspace(-1,0,17).';
            weights = [0.5;ones(15,1);0.5]/16;
            signed = IMBasisCollection({basis}).projectionOnGrid(z,weights,variable="u").metricMatrix;
            endpoint = evp.innerProduct().surfaceWeights(1);
            testCase.verifyEqual(signed(end,end),weights(end)+endpoint.coefficient*endpoint.c^2,AbsTol=1e-14)
            testCase.verifyError(@() IMBasisCollection({basis}).projectionOnGrid(z(1:end-1),weights(1:end-1),variable="u"),"IMProjection:MissingEndpointSample")
            evp = IMEigenvalueProblem(zDomain=[-1 0],p=1,q=0,r=1,surfaceBoundary=IMBoundaryCondition(a=1,b=0,c=0,d=1),bottomBoundary=IMBoundaryCondition.dirichlet());
            basis = solver.solveEVP(evp,nModes=3);
            testCase.verifyError(@() IMBasisCollection({basis}).projectionOnGrid(z,weights,variable="u"),"IMProjection:UnsupportedDerivativeTrace")
        end
        function signedAlignedAndMDAKeepDistinctCapabilities(testCase)
            N2 = @(z) 1e-4*ones(size(z));
            evp = IMInternalModes.geostrophicAPVModes(rotationRate=0,N2=N2,zDomain=[-1000 0],g0=-0.04,gd=0.03,surfaceBoundary="freeSurface");
            solver = IMSolverSpectral(nEVP=64);
            numerical = solver.solveEVP(evp,nModes=4);
            analytical = IMConstantStratificationSolution(rotationRate=0,N0=0.01,zDomain=[-1000 0]).internalModes(evp,nModes=4);
            for basis = {numerical,analytical}
                recipe = basis{1}.projectionRecipe(variable="G");
                testCase.verifyTrue(recipe.isAvailable)
                testCase.verifyTrue(recipe.supportsQuadratic)
                testCase.verifyFalse(recipe.supportsLeakage)
                testCase.verifyEqual(recipe.targetGramMatrix,basis{1}.gramMatrix(variable="G"),AbsTol=1e-12)
                testCase.verifyEqual(recipe.majorantGramMatrix,basis{1}.majorantGramMatrix(variable="G"),AbsTol=1e-12)
            end
            evp = IMInternalModes.meanDensityAnomalyModes(rotationRate=0,N2=N2,zDomain=[-1000 0],g0=-0.04,gd=0.03);
            basis = solver.solveEVP(evp,nModes=4);
            recipe = basis.projectionRecipe(variable="G");
            testCase.verifyTrue(recipe.isAvailable)
            testCase.verifyFalse(recipe.supportsQuadratic)
            testCase.verifyFalse(basis.projectionRecipe(variable="F").isAvailable)
        end
        function zeroAPVRecipesPreserveEndpointLabelsAndRejectScalarProjection(testCase)
            options = namedargs2cell(struct("N2",@(z) 1e-4*ones(size(z)),"zDomain",[-1000 0],"rotationRate",1e-4/(2*sind(33)),"latitude",33,"k",[1e-4 2e-4],"endpoints",["surface","bottom"]));
            problem = IMGeostrophicZeroAPVModes.atWavenumber(options{:});
            numerical = IMSolverSpectral(nEVP=48).solveGeostrophicZeroAPVModes(problem);
            analytical = IMConstantStratificationSolution(N0=0.01,zDomain=[-1000 0],rotationRate=1e-4/(2*sind(33)),latitude=33).geostrophicZeroAPVModesAtWavenumber([1e-4 2e-4],endpoints=["surface","bottom"]);
            for basis = {numerical,analytical}
                recipe = basis{1}.projectionRecipe(variable="F");
                testCase.verifyFalse(recipe.isAvailable)
                testCase.verifyEqual(recipe.columnLabels,["surface","bottom"])
                testCase.verifyEqual(recipe.provenance.requestedKappa,[1e-4 2e-4])
                testCase.verifyError(@() IMBasisCollection(basis).projectionOnGrid([-1000;0],[500;500],variable="F"),"IMScalarProjection:UnavailableProjection")
            end
        end
    end
end
