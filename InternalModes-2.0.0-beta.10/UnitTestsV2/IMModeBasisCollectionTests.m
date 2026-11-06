classdef IMModeBasisCollectionTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addProvider(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
        end
    end

    methods (Test)
        function wavePagesPreserveTheStandardConventionAndScientificIdentity(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            solver = IMSolverAnalytical();
            first = solver.solveEVP(IMInternalModes.waveModesAtWavenumber(background,k=1e-3),nModes=4).withNormalization("waveEnergy");
            second = solver.solveEVP(IMInternalModes.waveModesAtWavenumber(background,k=2e-3),nModes=4).withNormalization("waveEnergy");
            collection = IMBasisCollection({first,second},basisIndex=[2 1 2]);
            testCase.verifyEqual(collection.kappa,[2e-3 1e-3 2e-3]);
            testCase.verifyEqual(collection.metadata(1).family,first.family);
            testCase.verifyEqual(collection.metadata(1).columnLabels,first.columnLabels);
            testCase.verifyEqual(collection.metadata(1).normalization,"waveEnergy");
            testCase.verifyEqual(collection.metadata(1).variables,["F", "G"]);
            testCase.verifyEqual(collection.metadata(1).derivativeOrders,{[0 1],0:2});
            z = linspace(-1000,0,17).';
            pages = [3 2 3 1];
            columns = [3 1 3];
            F = collection.evaluate(z,variable="F",pages=pages,columns=columns,sampleChunkSize=4,pageChunkSize=2);
            Gz = collection.evaluate(z,variable="G",derivativeOrder=1,pages=pages,columns=columns,sampleChunkSize=3,pageChunkSize=1);
            testCase.verifyEqual(F,background.D*Gz,AbsTol=2e-14);
            for j = 1:numel(pages)
                basis = collection.bases{collection.basisIndex(pages(j))};
                expected = basis.F(z);
                testCase.verifyEqual(F(:,:,j),expected(:,columns),AbsTol=0);
            end
            derivativeF = collection.evaluate(z,variable="F",derivativeOrder=1);
            testCase.verifyEqual(derivativeF(:,:,1),second.evaluate(z,variable="F",derivativeOrder=1),AbsTol=1e-13);
        end

        function streamedNativePagesKeepRepeatedPositionsAndBoundedChunks(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            basis = IMSolverAnalytical().solveEVP(IMInternalModes.waveModesAtWavenumber(background,k=1e-3),nModes=3);
            collection = IMBasisCollection({basis},basisIndex=[1 1 1]);
            z = linspace(-1000,0,19).';
            pages = [3 1 3 2 1];
            columns = [3 1];
            output = zeros(numel(z),numel(columns),numel(pages));
            coverage = zeros(numel(z),numel(pages));
            callCount = 0;
            collection.evaluateChunks(z,@consume,variable="G",derivativeOrder=2,pages=pages,columns=columns,sampleChunkSize=4,pageChunkSize=2);
            expected = basis.Gzz(z);
            testCase.verifyEqual(output,repmat(expected(:,columns),1,1,numel(pages)),AbsTol=0);
            testCase.verifyEqual(coverage,ones(size(coverage)));
            testCase.verifyEqual(callCount,15);

            function consume(values,rows,positions,pageIDs)
                testCase.verifyLessThanOrEqual(numel(rows),4);
                testCase.verifyLessThanOrEqual(numel(positions),2);
                testCase.verifyEqual(pageIDs,pages(positions));
                output(rows,:,positions) = values;
                coverage(rows,positions) = coverage(rows,positions)+1;
                callCount = callCount+1;
            end
        end

        function endpointBasesUseNativeEvaluatorsAcrossCollectionPages(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            problem = IMGeostrophicZeroAPVModes.atWavenumber(background,1e-3,endpoints=["surface" "bottom"]);
            basis = IMSolverAnalytical().solveGeostrophicZeroAPVModes(problem);
            collection = IMBasisCollection({basis},basisIndex=[1 1]);
            testCase.verifyEqual(collection.metadata.columnKind,"endpoint");
            testCase.verifyEqual(collection.metadata.endpointLabels,["surface", "bottom"]);
            z = [-1000;-700;-100;0];
            for variable = ["F", "G"]
                expected = basis.evaluate(z,variable=variable);
                sampled = collection.evaluate(z,variable=variable,pages=[2 1 2],columns=[2 1],sampleChunkSize=2,pageChunkSize=2);
                testCase.verifyEqual(sampled,repmat(expected(:,[2 1]),1,1,3),AbsTol=0);
            end
            testCase.verifyError(@() IMBasisCollection({basis},sourcePage=2),"IMBasisCollection:InvalidSourcePage");
        end

        function emptyEndpointInventoryRemainsAnEmptyContinuousBasis(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            problem = IMGeostrophicZeroAPVModes.atWavenumber(background,1e-3,endpoints=strings(1,0));
            basis = IMSolverSpectral(nEVP=32).solveGeostrophicZeroAPVModes(problem);
            collection = IMBasisCollection({basis},basisIndex=[1 1]);
            testCase.verifyEqual(basis.nColumns,0);
            testCase.verifyEmpty(collection.metadata.columnLabels);
            testCase.verifySize(collection.evaluate([-1000;-500;0],variable="F"),[3 0 2]);
        end

        function neutralAPVProjectionPreservesNativeNormalization(testCase)
            background = IMBackground(N2=@(z) 1e-4*double(z < -200),zDomain=[-1000 0],surfaceCondition="freeSurface");
            depth = IMNeutralLabelModesBasis(background,neutralIntervals=[-200 0],g0=0,gd=Inf,nModesPerInterval=3);
            maximum = IMNeutralLabelModesBasis(background,neutralIntervals=[-200 0],g0=0,gd=Inf,nModesPerInterval=3,normalization="maxAbsG");
            collection = IMBasisCollection({depth,maximum},basisIndex=[2 1 2]);
            testCase.verifyEqual(collection.metadata(1).family,"neutralLabel");
            testCase.verifyEqual(collection.metadata(1).variables,["G", "APV"]);
            testCase.verifyEqual(collection.metadata(1).derivativeOrders,{0:2,0});
            z = linspace(-200,0,401).';
            weights = .5*ones(size(z));
            weights([1 end]) = .25;
            coefficients = [2;-.5;1.25];
            for page = [1 2]
                basis = collection.bases{collection.basisIndex(page)};
                projection = collection.projectionOnGrid(z,weights,variable="APV",page=page);
                explicit = basis.projectionOnGrid(z,weights,variable="APV");
                testCase.verifyEqual(projection.physicalDualMatrix,explicit.physicalDualMatrix,AbsTol=0);
                testCase.verifyEqual(projection.columnLabels,basis.columnLabels);
                testCase.verifyEqual(projection.physicalDualMatrix*basis.APV(z)*coefficients,coefficients,AbsTol=2e-13);
                selected = collection.projectionOnGrid(z,weights,variable="APV",page=page,columns=[3 1]);
                samples = basis.APV(z);
                testCase.verifyEqual(selected.physicalDualMatrix*samples(:,[3 1]),eye(2),AbsTol=2e-13);
            end
            sampled = collection.evaluate(z,variable="APV",pages=[3 1 2],columns=[3 1]);
            testCase.verifyEqual(sampled(:,:,1),sampled(:,:,2),AbsTol=0);
            testCase.verifyError(@() collection.projectionOnGrid(z,weights,variable="G"),"IMScalarProjection:UnavailableProjection");
        end

        function emptyNeutralViewsAndScaledSelectionsPreserveDualMetadata(testCase)
            background = IMBackground(N2=@(z) 1e-4*double(z < -200),zDomain=[-1000 0],surfaceCondition="freeSurface");
            maximum = IMNeutralLabelModesBasis(background,neutralIntervals=[-200 0],nModesPerInterval=3,normalization="maxAbsG");
            requestedEmpty = IMNeutralLabelModesBasis(background,neutralIntervals=[-200 0],nModesPerInterval=0);
            selectedEmpty = maximum.selectColumns([]);
            collection = IMBasisCollection({requestedEmpty,selectedEmpty},basisIndex=[2 1 2]);
            z = linspace(-200,0,401).';
            weights = .5*ones(size(z)); weights([1 end]) = .25;
            for basis = {requestedEmpty,selectedEmpty}
                projection = basis{1}.projectionOnGrid(z,weights,variable="APV");
                recipe = basis{1}.projectionRecipe();
                testCase.verifySize(projection.physicalDualMatrix,[0 numel(z)]);
                testCase.verifyEqual(projection.majorantGramMatrix,zeros(0));
                testCase.verifyEqual(projection.provenance.continuousAPVGramMatrix,zeros(0));
                testCase.verifyEqual(recipe.provenance.continuousAPVGramMatrix,zeros(0));
                testCase.verifyEmpty(projection.columnLabels);
                testCase.verifySize(projection.projectPhysical(ones(numel(z),2)),[0 2]);
            end
            for derivativeOrder = 0:2
                testCase.verifySize(collection.evaluate(z,variable="G",derivativeOrder=derivativeOrder,pages=[3 1 3 2],sampleChunkSize=100,pageChunkSize=2),[numel(z) 0 4]);
            end
            testCase.verifySize(collection.evaluate(z,variable="APV",pages=[3 1 3 2]),[numel(z) 0 4]);
            for page = [3 1 3 2]
                projection = collection.projectionOnGrid(z,weights,variable="APV",page=page);
                testCase.verifySize(projection.physicalDualMatrix,[0 numel(z)]);
                testCase.verifyEqual(projection.provenance.continuousAPVGramMatrix,zeros(0));
            end
            selected = maximum.selectColumns([3 1]);
            expected = diag(background.f0^2*100*(pi*[3 1]/200).^2);
            projection = selected.projectionOnGrid(z,weights,variable="APV");
            recipe = selected.projectionRecipe();
            testCase.verifyEqual(projection.majorantGramMatrix,expected,RelTol=5e-15);
            testCase.verifyEqual(projection.provenance.continuousAPVGramMatrix,expected,RelTol=5e-15);
            testCase.verifyEqual(recipe.majorantGramMatrix,expected,RelTol=5e-15);
            testCase.verifyEqual(recipe.provenance.continuousAPVGramMatrix,expected,RelTol=5e-15);
            testCase.verifyEqual(projection.columnLabels,maximum.columnLabels([3 1]));
            testCase.verifyEqual(projection.physicalDualMatrix*selected.APV(z),eye(2),AbsTol=2e-13);
        end

        function emptyNeutralViewDoesNotPrepareDiscardedAPVNorms(testCase)
            background = IMBackground(N2=@(z) zeros(size(z)),zDomain=[-1000 0],rotationRate=1e-170,surfaceCondition="freeSurface");
            full = IMNeutralLabelModesBasis(background,neutralIntervals=[-1000 0],nModesPerInterval=3);
            selected = full.selectColumns([]);
            requested = IMNeutralLabelModesBasis(background,neutralIntervals=[-1000 0],nModesPerInterval=0);
            z = [-900;-500;-100];
            weights = [300;400;300];
            testCase.verifyError(@() full.projectionOnGrid(z,weights,variable="APV"),"IMNeutralLabelModesBasis:UnrepresentableProjection");
            collection = IMBasisCollection({selected,requested},basisIndex=[2 1 2]);
            for page = 1:3
                projection = collection.projectionOnGrid(z,weights,variable="APV",page=page);
                testCase.verifySize(projection.physicalDualMatrix,[0 numel(z)]);
                testCase.verifyEqual(projection.physicalDualStatus,"available");
                testCase.verifySize(projection.projectPhysical(ones(numel(z),2)),[0 2]);
                testCase.verifyEqual(projection.provenance.continuousAPVGramMatrix,zeros(0));
                testCase.verifyEmpty(projection.recovery);
            end
            testCase.verifyError(@() selected.projectionOnGrid(z,zeros(size(z)),variable="APV"),"IMNeutralLabelModesBasis:InvalidWeights");
            testCase.verifyError(@() selected.projectionOnGrid([z(1);z(1)],weights(1:2),variable="APV"),"IMNeutralLabelModesBasis:InvalidGrid");
        end

        function collectionPreservesAuxiliaryWaveAssessmentCapabilities(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            basis = IMSolverAnalytical().solveEVP(IMInternalModes.waveModesAtWavenumber(background,k=1e-3),nModes=3);
            collection = IMBasisCollection({basis});
            z = linspace(-1000,0,11).';
            weights = 100*ones(size(z));
            weights([1 end]) = 50;
            testCase.verifyError(@() collection.projectionOnGrid(z,weights,variable="G",prepareRecovery=true),"IMScalarProjection:UnavailableProjection");
            assessment = collection.assess(z,weights,variable="G");
            projection = assessment.projection;
            testCase.verifyEqual(projection.columnLabels,basis.columnLabels);
            testCase.verifyEqual(projection.targetGramMatrix,basis.gramMatrix(variable="G"),AbsTol=1e-12);
            testCase.verifyEqual(projection.recovery.transformForward(basis.G(z)),eye(3),AbsTol=1e-12);
        end
    end
end
