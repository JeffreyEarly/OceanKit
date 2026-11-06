classdef IMEvaluationConsistencyTests < matlab.unittest.TestCase
    % Verify the common scalar, collection and streamed evaluation contract.
    methods (TestClassSetup)
        function configurePath(testCase)
            original = path;
            root = fileparts(fileparts(mfilename("fullpath")));
            addpath(root,fullfile(root,"UnitTestsV2","TestSupport"));
            testCase.addTeardown(@() path(original));
        end
    end

    methods (Test)
        function explicitEmptyPagesDoNotEvaluateAnyBasis(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            collection = IMSolverAnalytical().solveWaveModesAtWavenumbers(background,[0 1e-4],nModes=[2 3]);
            z = [-900;-500;-100];
            testCase.verifySize(collection.evaluate(z,variable="F",pages=[]),[3 0 0]);
            testCase.verifySize(collection.evaluate(z,variable="F",pages=[],columns=[3 1 3]),[3 3 0]);
            testCase.verifySize(collection.evaluate(z,variable="F",pages=[],columns=[]),[3 0 0]);
            callbacks = 0;
            collection.evaluateChunks(z,@consume,variable="F",pages=[]);
            testCase.verifyEqual(callbacks,0);
            testCase.verifyError(@() collection.evaluate(z,variable="F"),"IMBasisCollection:HeterogeneousColumns");
            testCase.verifyError(@() collection.evaluate(z,variable="F",pages=[],columns=[0 1]),"IMBasisCollection:InvalidColumn");
            function consume(varargin)
                callbacks = callbacks+1;
            end
        end

        function variableIsExplicitForEveryFamily(testCase)
            background = IMBackground.constantStratification(surfaceCondition="freeSurface",N0=.01,zDomain=[-1000 0]);
            problem = IMInternalModes.waveModesAtWavenumber(background,k=2e-4);
            numeric = IMSolverSpectral(nEVP=32).solveEVP(problem,nModes=2);
            exact = IMSolverAnalytical().solveEVP(problem,nModes=2);
            neutral = IMBackground(surfaceCondition="freeSurface",N2=@(z) zeros(size(z)),zDomain=background.zDomain);
            labels = IMNeutralLabelModesBasis(neutral,neutralIntervals=background.zDomain,nModesPerInterval=2);
            critical = IMBackground(surfaceCondition="freeSurface",N2=@(z) background.f0^2*ones(size(z)),zDomain=background.zDomain);
            zeroDepth = IMZeroDepthWaveModesBasis(critical,2e-4,criticalIntervals=background.zDomain,nModesPerInterval=2);
            endpointProblem = IMGeostrophicZeroAPVModes.atWavenumber(background,2e-4,endpoints="surface");
            endpoint = IMSolverSpectral(nEVP=32).solveGeostrophicZeroAPVModes(endpointProblem);
            exactEndpoint = IMSolverAnalytical().solveGeostrophicZeroAPVModes(endpointProblem);
            callbacks = 0;
            for item = {numeric,exact,labels,zeroDepth,endpoint,exactEndpoint}
                basis = item{1};
                testCase.verifyError(@() basis.evaluate(-500),string(class(basis))+":MissingVariable");
                collection = IMBasisCollection({basis});
                testCase.verifyError(@() collection.evaluate(-500),"IMBasisCollection:MissingVariable");
                testCase.verifyError(@() collection.evaluateChunks(-500,@consume),"IMBasisCollection:MissingVariable");
                for alias = ["Fz","Gz","Gzz"]
                    testCase.verifyError(@() basis.evaluate(-500,variable=alias),string(class(basis))+":UnsupportedVariable");
                    testCase.verifyError(@() collection.evaluate(-500,variable=alias),"IMBasisCollection:UnsupportedVariable");
                end
            end
            testCase.verifyEqual(callbacks,0);
            function consume(varargin)
                callbacks = callbacks+1;
            end
        end

        function normalizedCollectionsAgreeWithScalarBackends(testCase)
            background = IMBackground.constantStratification(surfaceCondition="freeSurface",N0=.01,zDomain=[-1000 0]);
            problem = IMInternalModes.waveModesAtWavenumber(background,k=2e-4);
            solvers = {IMSolverSpectral(nEVP=40),IMSolverFiniteDifference(z=linspace(-1000,0,65).'),IMSolverAnalytical()};
            bases = cell(size(solvers));
            for i = 1:numel(solvers)
                bases{i} = solvers{i}.solveEVP(problem,nModes=3).addNormalization("custom",@(~,j) j+2);
            end
            collection = IMBasisCollection(bases,basisIndex=[3 1 2 1]);
            z = linspace(-1000,0,19).';
            columns = [3 1 3];
            for variable = ["F","G"]
                original = collection.evaluate(z,variable=variable);
                for order = 0:1
                    values = collection.evaluate(z,variable=variable,derivativeOrder=order,normalization="custom",columns=columns,sampleChunkSize=4,pageChunkSize=1);
                    for page = 1:numel(collection.basisIndex)
                        basis = bases{collection.basisIndex(page)};
                        expected = basis.evaluate(z,variable=variable,derivativeOrder=order,normalization="custom");
                        testCase.verifyEqual(values(:,:,page),expected(:,columns),AbsTol=2e-12,RelTol=2e-12);
                    end
                end
                testCase.verifyEqual(collection.evaluate(z,variable=variable),original,AbsTol=0);
            end
            testCase.verifyEqual(string({collection.metadata.normalization}),repmat("waveEnergy",1,3));
            testCase.verifyEqual(collection.bases,bases);
        end

        function streamedOverridesPreserveSelectionsAndDerivatives(testCase)
            background = IMBackground.constantStratification(surfaceCondition="freeSurface",N0=.01,zDomain=[-1000 0]);
            problem = IMInternalModes.waveModesAtWavenumber(background,k=2e-4);
            spectral = IMSolverSpectral(nEVP=32).solveEVP(problem,nModes=3);
            exact = IMSolverAnalytical().solveEVP(problem,nModes=3);
            collection = IMBasisCollection({spectral,exact},basisIndex=[2 1 2]);
            z = linspace(-1000,0,17).';
            pages = [3 2 3 1]; columns = [2 1 2];
            actual = zeros(numel(z),numel(columns),numel(pages));
            coverage = zeros(numel(z),numel(pages));
            collection.evaluateChunks(z,@consume,variable="G",derivativeOrder=2,normalization="waveEnergy",pages=pages,columns=columns,sampleChunkSize=5,pageChunkSize=2);
            expected = collection.evaluate(z,variable="G",derivativeOrder=2,normalization="waveEnergy",pages=pages,columns=columns);
            testCase.verifyEqual(actual,expected,AbsTol=1e-12,RelTol=1e-12);
            testCase.verifyEqual(coverage,ones(size(coverage)));
            function consume(values,rows,positions,pageIDs)
                testCase.verifyEqual(pageIDs,pages(positions));
                actual(rows,:,positions) = values;
                coverage(rows,positions) = coverage(rows,positions)+1;
            end
        end

        function nonlinearCompanionIsRecoveredBeforeNormalization(testCase)
            problem = IMInternalModes(rotationRate=0,name="customDiagnostic",formulation="G",N2=@(z) ones(size(z)),zDomain=[-1 0],p=1,q=0,r=1,FfromGz=@(~,dGdz,~,~) dGdz.^2+1);
            basis = IMSolverSpectral(nEVP=32).solveEVP(problem,nModes=3).addNormalization("triple",@(~,~) 3);
            collection = IMBasisCollection({basis},basisIndex=[1 1]);
            z = linspace(-1,0,13).';
            expected = basis.F(z,normalization="triple");
            actual = collection.evaluate(z,variable="F",normalization="triple",sampleChunkSize=4);
            testCase.verifyEqual(actual,repmat(expected,1,1,2),AbsTol=1e-11,RelTol=1e-12);
        end

        function genericScalarAndSpecialFamiliesAcceptOverrides(testCase)
            solver = IMSolverSpectral(nEVP=32);
            problem = IMEigenvalueProblem(zDomain=[-1 0],p=1,q=0,r=1,surfaceBoundary=IMBoundaryCondition.dirichlet(),bottomBoundary=IMBoundaryCondition.dirichlet());
            scalar = solver.solveEVP(problem,nModes=3).addNormalization("custom",@(~,j) j+1);
            z = linspace(-1,0,13).';
            collection = IMBasisCollection({scalar});
            testCase.verifyError(@() scalar.evaluate(z),"IMBasisSet:MissingVariable");
            for order = 0:1
                testCase.verifyEqual(collection.evaluate(z,variable="u",derivativeOrder=order,normalization="custom"),scalar.evaluate(z,variable="u",derivativeOrder=order,normalization="custom"),AbsTol=1e-12);
            end
            background = IMBackground(surfaceCondition="freeSurface",N2=@(z) zeros(size(z)),zDomain=[-1 0]);
            labels = IMNeutralLabelModesBasis(background,neutralIntervals=[-1 0],nModesPerInterval=3).addNormalization("custom",@(~,j) j+1);
            critical = IMBackground(surfaceCondition="freeSurface",N2=@(z) background.f0^2*ones(size(z)),zDomain=[-1 0]);
            waves = IMZeroDepthWaveModesBasis(critical,2e-4,criticalIntervals=[-1 0],nModesPerInterval=3).addNormalization("custom",@(~,j) j+1);
            for entry = {{labels,"APV",0},{labels,"G",2},{waves,"F",1},{waves,"G",2}}
                item = entry{1}; basis = item{1}; variable = item{2}; order = item{3};
                collection = IMBasisCollection({basis},basisIndex=[1 1]);
                expected = basis.evaluate(z,variable=variable,derivativeOrder=order,normalization="custom");
                actual = collection.evaluate(z,variable=variable,derivativeOrder=order,normalization="custom",sampleChunkSize=4,pageChunkSize=1);
                testCase.verifyEqual(actual,repmat(expected,1,1,2),AbsTol=1e-11,RelTol=1e-12);
            end
        end

        function invalidNormalizationFailsBeforeAnyCallback(testCase)
            background = IMBackground.constantStratification(surfaceCondition="freeSurface",N0=.01,zDomain=[-1000 0]);
            problem = IMInternalModes.waveModesAtWavenumber(background,k=2e-4);
            basis = IMSolverSpectral(nEVP=32).solveEVP(problem,nModes=2);
            custom = basis.addNormalization("onlyFirst",@(~,~) 2);
            collection = IMBasisCollection({custom,basis});
            callbacks = 0;
            testCase.verifyError(@() collection.evaluateChunks([-900;-500],@consume,variable="G",normalization="onlyFirst",sampleChunkSize=1,pageChunkSize=1),"IMBasisSet:UnsupportedNormalization");
            testCase.verifyEqual(callbacks,0);
            function consume(varargin)
                callbacks = callbacks+1;
            end
        end

        function metadataUsesActualDerivativeCapabilities(testCase)
            background = IMBackground.constantStratification(surfaceCondition="freeSurface",N0=1,zDomain=[-1 0],rotationRate=0);
            ordinary = IMSolverSpectral(nEVP=32).solveEVP(IMInternalModes.waveModesAtWavenumber(background,k=2e-4),nModes=2);
            problem = IMInternalModes(rotationRate=0,name="customDiagnostic",formulation="G",N2=@(z) ones(size(z)),zDomain=[-1 0],p=1,q=0,r=1,FfromGz=@(~,dGdz,~,~) dGdz.^2+1);
            custom = IMSolverSpectral(nEVP=32).solveEVP(problem,nModes=2);
            collection = IMBasisCollection({ordinary,custom});
            testCase.verifyEqual(collection.metadata(2).derivativeOrders{1},0);
            callbacks = 0;
            testCase.verifyError(@() collection.evaluateChunks([-.9;-.5],@consume,variable="F",derivativeOrder=1,sampleChunkSize=1,pageChunkSize=1),"IMBasisCollection:UnsupportedDerivative");
            testCase.verifyEqual(callbacks,0);
            exact = IMAnalyticalInternalModesBasis(evp=problem,h=1/pi^2,modeNumber=1,N2=problem.N2,rawVariableFunction=@(~,variable,z) (variable=="G").*sin(pi*(z+1))+(variable=="F").*(pi^2*cos(pi*(z+1)).^2+1),rawUzFunction=@(~,z) pi*cos(pi*(z+1)));
            exactCollection = IMBasisCollection({exact});
            testCase.verifyEqual(exactCollection.metadata.derivativeOrders,{0,0:1});
            testCase.verifyError(@() exactCollection.evaluateChunks([-.9;-.5],@consume,variable="G",derivativeOrder=2),"IMBasisCollection:UnsupportedDerivative");
            testCase.verifyEqual(callbacks,0);
            function consume(varargin)
                callbacks = callbacks+1;
            end
        end

        function endpointNormalizationPreservesCoordinateMeaning(testCase)
            background = IMBackground.constantStratification(surfaceCondition="freeSurface",N0=.01,zDomain=[-1000 0]);
            problem = IMGeostrophicZeroAPVModes.atWavenumber(background,[1e-4 2e-4]);
            numerical = IMSolverSpectral(nEVP=40).solveGeostrophicZeroAPVModes(problem);
            exact = IMSolverAnalytical().solveGeostrophicZeroAPVModes(problem);
            z = linspace(-1000,0,13).';
            callbacks = 0;
            for entry = {numerical,exact,numerical.orthogonalizeGeneralizedEnergy(g0=0,gd=0),exact.orthogonalizeGeneralizedEnergy(g0=0,gd=0)}
                basis = entry{1};
                collection = IMBasisCollection({basis},basisIndex=[1 1 1],sourcePage=[2 1 2]);
                actual = collection.evaluate(z,variable="G",derivativeOrder=1,normalization="rmsF",sampleChunkSize=4);
                expected = basis.evaluate(z,variable="G",derivativeOrder=1,normalization="rmsF",pages=[2 1 2]);
                testCase.verifyEqual(actual,expected,AbsTol=1e-12,RelTol=1e-12);
                testCase.verifyError(@() collection.evaluateChunks(z,@consume,variable="F",normalization="notARule",sampleChunkSize=1),string(class(basis))+":UnsupportedNormalization");
                testCase.verifyError(@() basis.evaluate(z,variable="F",normalization="notARule"),string(class(basis))+":UnsupportedNormalization");
            end
            testCase.verifyEqual(callbacks,0);
            function consume(varargin)
                callbacks = callbacks+1;
            end
        end
    end
end
