classdef IMProjectionOnGridTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function configurePath(testCase)
            originalPath = path;
            addpath(fileparts(fileparts(mfilename("fullpath"))));
            testCase.addTeardown(@() path(originalPath));
        end
    end
    methods (Test)
        function scalarAndAlignedRecipesPreserveExistingTransforms(testCase)
            solver = IMSolverSpectral(nEVP=48);
            scalarEVP = IMEigenvalueProblem(zDomain=[-1 0],p=1,q=0,r=1,surfaceBoundary=IMBoundaryCondition.dirichlet(),bottomBoundary=IMBoundaryCondition.dirichlet());
            scalar = solver.solveEVP(scalarEVP,nModes=4);
            z = linspace(-1,0,25).';
            weights = [0.5;ones(23,1);0.5]/24;
            old = scalar.discreteTransform(z=z,weights=weights,nModes=3);
            collection = IMBasisCollection({scalar});
            projection = collection.projectionOnGrid(z,weights,columns=1:3);
            testCase.verifyEqual(projection.sampledBasis,old.inverseMatrix,AbsTol=1e-12)
            testCase.verifyEqual(projection.metricMatrix,old.metricMatrix,AbsTol=1e-12)
            testCase.verifyEqual(projection.targetGramMatrix,old.targetGramMatrix,AbsTol=1e-12)
            testCase.verifyEqual(projection.forwardMatrix,old.forwardMatrix,AbsTol=1e-11)
            evp = IMInternalModes.hydrostaticFModes(N2=@(z) ones(size(z)),zDomain=[-1 0]);
            aligned = solver.solveEVP(evp,nModes=4);
            old = aligned.discreteTransform(z=z,weights=weights,nModes=4,gramTolerance=1);
            for variable = ["F","G"]
                projection = IMBasisCollection({aligned}).projectionOnGrid(z,weights,variable=variable);
                testCase.verifyEqual(projection.sampledBasis,old.inverseMatrix(variable=variable),AbsTol=1e-12)
                testCase.verifyEqual(projection.forwardMatrix,old.forwardMatrix(variable=variable),AbsTol=1e-10)
                testCase.verifyEqual(projection.activeColumnMask,old.activeModeMask(variable=variable))
            end
        end
        function analyticalRecipesRetainNormalizationAndReferenceProvenance(testCase)
            evp = IMInternalModes.hydrostaticGModes(N2=@(z) 1e-4*ones(size(z)),zDomain=[-1000 0]);
            basis = IMConstantStratificationSolution(N0=0.01,zDomain=[-1000 0]).internalModes(evp,nModes=4);
            collection = IMBasisCollection({basis});
            z = linspace(-1000,0,21).';
            weights = [0.5;ones(19,1);0.5]*50;
            projection = collection.projectionOnGrid(z,weights,variable="G");
            expected = basis.G(z);
            basis.normalization = "wMax";
            testCase.verifyEqual(collection.evaluate(z,variable="G"),expected,AbsTol=1e-12)
            testCase.verifyEqual(collection.evaluate(z,variable="G",columns=[3 1]),expected(:,[3 1]),AbsTol=1e-12)
            testCase.verifyEqual(projection.provenance.referencePointCount,1024)
            testCase.verifyEqual(projection.provenance.referenceConvergence,"unverified")
            testCase.verifyEqual(projection.provenance.solveAccuracy,"unverified")
            testCase.verifyEqual(projection.prefix(2).provenance,projection.provenance)
            testCase.verifyEqual(collection.assess(z,weights,variable="G").identity.projectionProvenance,projection.provenance)
            testCase.verifyError(@() collection.projectionOnGrid(z,weights,variable="G",columns=[1 1]),"IMProjection:InvalidColumns")
        end
        function signedEndpointMetricsAndMissingTracesAreExplicit(testCase)
            solver = IMSolverSpectral(nEVP=48);
            evp = IMEigenvalueProblem(zDomain=[-1 0],p=1,q=0,r=1,surfaceBoundary=IMBoundaryCondition(a=0,b=1,c=1,d=0),bottomBoundary=IMBoundaryCondition.dirichlet());
            basis = solver.solveEVP(evp,nModes=3);
            z = linspace(-1,0,17).';
            weights = [0.5;ones(15,1);0.5]/16;
            signed = IMBasisCollection({basis}).projectionOnGrid(z,weights).metricMatrix;
            endpoint = evp.innerProduct().surfaceWeights(1);
            testCase.verifyEqual(signed(end,end),weights(end)+endpoint.coefficient*endpoint.c^2,AbsTol=1e-14)
            testCase.verifyError(@() IMBasisCollection({basis}).projectionOnGrid(z(1:end-1),weights(1:end-1)),"IMProjection:MissingEndpointSample")
            evp = IMEigenvalueProblem(zDomain=[-1 0],p=1,q=0,r=1,surfaceBoundary=IMBoundaryCondition(a=1,b=0,c=0,d=1),bottomBoundary=IMBoundaryCondition.dirichlet());
            basis = solver.solveEVP(evp,nModes=3);
            testCase.verifyError(@() IMBasisCollection({basis}).projectionOnGrid(z,weights),"IMProjection:UnsupportedDerivativeTrace")
        end
        function signedAlignedAndMDAKeepDistinctCapabilities(testCase)
            N2 = @(z) 1e-4*ones(size(z));
            evp = IMInternalModes.geostrophicAPVModes(N2=N2,zDomain=[-1000 0],g0=-0.04,gd=0.03);
            solver = IMSolverSpectral(nEVP=64);
            numerical = solver.solveEVP(evp,nModes=4);
            analytical = IMConstantStratificationSolution(N0=0.01,zDomain=[-1000 0]).internalModes(evp,nModes=4);
            for basis = {numerical,analytical}
                recipe = basis{1}.projectionRecipe(variable="G");
                testCase.verifyTrue(recipe.isAvailable)
                testCase.verifyTrue(recipe.supportsQuadratic)
                testCase.verifyFalse(recipe.supportsLeakage)
                testCase.verifyEqual(recipe.targetGramMatrix,basis{1}.gramMatrix(variable="G"),AbsTol=1e-12)
                testCase.verifyEqual(recipe.majorantGramMatrix,basis{1}.majorantGramMatrix(variable="G"),AbsTol=1e-12)
            end
            evp = IMInternalModes.meanDensityAnomalyModes(N2=N2,zDomain=[-1000 0],g0=-0.04,gd=0.03);
            basis = solver.solveEVP(evp,nModes=4);
            recipe = basis.projectionRecipe(variable="G");
            testCase.verifyTrue(recipe.isAvailable)
            testCase.verifyFalse(recipe.supportsQuadratic)
            testCase.verifyFalse(basis.projectionRecipe(variable="F").isAvailable)
        end
        function zeroAPVRecipesPreserveEndpointLabelsAndRejectScalarProjection(testCase)
            options = namedargs2cell(struct("N2",@(z) 1e-4*ones(size(z)),"zDomain",[-1000 0],"f0",1e-4,"k",[1e-4 2e-4],"endpoints",["surface","bottom"]));
            problem = IMGeostrophicZeroAPVModes.atWavenumber(options{:});
            numerical = IMSolverSpectral(nEVP=48).solveGeostrophicZeroAPVModes(problem);
            analytical = IMConstantStratificationSolution(N0=0.01,zDomain=[-1000 0],f0=1e-4).geostrophicZeroAPVModesAtWavenumber([1e-4 2e-4],endpoints=["surface","bottom"]);
            for basis = {numerical,analytical}
                recipe = basis{1}.projectionRecipe(variable="F");
                testCase.verifyFalse(recipe.isAvailable)
                testCase.verifyEqual(recipe.columnLabels,["surface","bottom"])
                testCase.verifyEqual(recipe.provenance.requestedKappa,[1e-4 2e-4])
                testCase.verifyError(@() IMBasisCollection(basis).projectionOnGrid([-1000;0],[500;500],variable="F"),"IMProjection:UnavailableProjection")
            end
        end
    end
end
