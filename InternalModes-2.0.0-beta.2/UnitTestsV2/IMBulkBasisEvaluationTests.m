classdef IMBulkBasisEvaluationTests < matlab.unittest.TestCase
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
        function restorePath(testCase)
            path(testCase.originalPath);
        end
    end
    methods (Test)
        function spectralPagesAndDerivativesAgreeWithScalarEvaluation(testCase)
            collection = testCase.waveCollection();
            z = linspace(-1000,0,31).';
            requested = [3 1 2 1];
            columns = [3 1 3];
            for variable = ["F","G"]
                values = collection.evaluate(z,variable=variable,pages=requested,columns=columns,sampleChunkSize=7,pageChunkSize=3);
                for i = 1:numel(requested)
                    basis = collection.bases{collection.basisIndex(requested(i))};
                    expected = basis.(variable)(z);
                    testCase.verifyEqual(values(:,:,i),expected(:,columns),RelTol=2e-12,AbsTol=2e-12);
                end
            end
            derivatives = collection.evaluate(z,variable="G",derivativeOrder=1,sampleChunkSize=6,pageChunkSize=2);
            for i = 1:numel(collection.basisIndex)
                expected = collection.bases{collection.basisIndex(i)}.uz(z);
                testCase.verifyEqual(derivatives(:,:,i),expected,RelTol=2e-12,AbsTol=2e-12);
            end
        end
        function streamingHasBoundedChunksAndUnambiguousRepeatedPositions(testCase)
            collection = testCase.waveCollection();
            z = linspace(-1000,0,23).';
            requested = [3 1 3 2 1];
            output = zeros(23,3,5);
            coverage = zeros(23,5);
            callCount = 0;
            collection.evaluateChunks(z,@consume,variable="F",pages=requested,sampleChunkSize=5,pageChunkSize=2);
            expected = collection.evaluate(z,variable="F",pages=requested,sampleChunkSize=100,pageChunkSize=100);
            testCase.verifyEqual(output,expected,RelTol=2e-12,AbsTol=2e-12);
            testCase.verifyEqual(coverage,ones(size(coverage)));
            testCase.verifyEqual(callCount,15);
            function consume(values,rows,positions,pageIDs)
                testCase.verifyLessThanOrEqual(numel(rows),5);
                testCase.verifyLessThanOrEqual(numel(positions),2);
                testCase.verifyEqual(pageIDs,requested(positions));
                output(rows,:,positions) = values;
                coverage(rows,positions) = coverage(rows,positions)+1;
                callCount = callCount+1;
            end
        end
        function scalarAndAnalyticalFallbackPreserveNormalization(testCase)
            solver = IMSolverSpectral(nEVP=40);
            scalarEVP = IMEigenvalueProblem(zDomain=[-1 0],p=1,q=0,r=1,surfaceBoundary=IMBoundaryCondition.dirichlet(),bottomBoundary=IMBoundaryCondition.dirichlet());
            scalar = solver.solveEVP(scalarEVP,nModes=3);
            collection = IMBasisCollection({scalar});
            z = linspace(-1,0,11).';
            testCase.verifyEqual(collection.evaluate(z,sampleChunkSize=3),scalar.u(z),AbsTol=1e-13);
            testCase.verifyEqual(collection.evaluate(z,derivativeOrder=1,sampleChunkSize=3),scalar.uz(z),AbsTol=1e-12);
            N2 = @(z) 2.5e-5*ones(size(z));
            evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=[-1000 0]);
            solution = IMConstantStratificationSolution(N0=5e-3,zDomain=[-1000 0]);
            analytical = solution.internalModes(evp,nModes=3);
            mean = solver.solveEVP(IMInternalModes.meanDensityAnomalyModes(N2=N2,zDomain=[-1000 0],g0=0,gd=0),nModes=3);
            for basis = {analytical,mean}
                collection = IMBasisCollection(basis,basisIndex=[1 1]);
                for variable = ["F","G"]
                    values = collection.evaluate(1000*z,variable=variable,sampleChunkSize=3,pageChunkSize=1);
                    testCase.verifyEqual(values(:,:,1),basis{1}.(variable)(1000*z),AbsTol=2e-12);
                    testCase.verifyEqual(values(:,:,1),values(:,:,2),AbsTol=0);
                end
            end
        end
        function boundarySelectionPreservesSourcePagesAndRotations(testCase)
            N2 = @(z) 2.5e-5*ones(size(z));
            k = [1e-4 2e-4 3e-4];
            problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=[-1000 0],f0=1e-4,k=k);
            numerical = IMSolverSpectral(nEVP=48).solveGeostrophicZeroAPVModes(problem);
            constant = IMConstantStratificationSolution(N0=5e-3,zDomain=[-1000 0],f0=1e-4);
            exponential = IMExponentialStratificationSolution(N0=5e-3,b=1000,zDomain=[-1000 0],f0=1e-4);
            z = linspace(-1000,0,17).';
            for entry = {numerical,constant.geostrophicZeroAPVModesAtWavenumber(k),exponential.geostrophicZeroAPVModesAtWavenumber(k)}
                basis = entry{1}.rotateBoundaryDepth();
                collection = IMBasisCollection({basis},basisIndex=[1 1 1],sourcePage=[3 1 2]);
                for variable = ["F","G"]
                    full = basis.(variable)(z);
                    selected = basis.(variable)(z,pages=[3 1 3]);
                    testCase.verifyEqual(selected,full(:,:,[3 1 3]),RelTol=3e-14,AbsTol=1e-15);
                    values = collection.evaluate(z,variable=variable,pages=[2 1 2],sampleChunkSize=4,pageChunkSize=2);
                    testCase.verifyEqual(values,full(:,:,[1 3 1]),RelTol=1e-12,AbsTol=2e-12);
                end
                derivative = collection.evaluate(z,variable="F",derivativeOrder=1,sampleChunkSize=4,pageChunkSize=2);
                expected = -(basis.N2(z)./basis.g).*basis.G(z,pages=[3 1 2]);
                testCase.verifyEqual(derivative,expected,RelTol=1e-12,AbsTol=2e-12);
            end
        end
        function customDiagnosticRecoveryPrecedesNormalization(testCase)
            evp = IMInternalModes(name="customDiagnostic",formulation="G",N2=@(z) ones(size(z)),zDomain=[-1 0],p=1,q=0,r=1,FfromGz=@(z,dGdz,h,ctx) dGdz.^2+1);
            basis = IMSolverSpectral(nEVP=32).solveEVP(evp,nModes=3);
            basis = basis.addNormalization("triple",@(b,j) 3*b.innerProductNormFactor(j));
            basis.normalization = "triple";
            collection = IMBasisCollection({basis},basisIndex=[1 1]);
            z = linspace(-1,0,17).';
            expected = basis.F(z);
            values = collection.evaluate(z,variable="F",sampleChunkSize=5,pageChunkSize=2);
            testCase.verifyEqual(values,repmat(expected,1,1,2),RelTol=1e-12,AbsTol=1e-11);
        end
        function invalidRequestsFailBeforeCallingConsumer(testCase)
            collection = testCase.waveCollection();
            callbacks = 0;
            testCase.verifyError(@() collection.evaluateChunks([-100;1],@consume,sampleChunkSize=1),"IMBasisCollection:OutsideDomain");
            testCase.verifyError(@() collection.evaluateChunks(-100,@consume,variable="F",derivativeOrder=1),"IMBasisCollection:UnsupportedDerivative");
            testCase.verifyEqual(callbacks,0);
            function consume(varargin)
                callbacks = callbacks+1;
            end
        end
    end
    methods (Static, Access = private)
        function collection = waveCollection()
            N2 = @(z) 2.5e-5*ones(size(z));
            solver = IMSolverSpectral(nEVP=48);
            bases = cell(1,2);
            for i = 1:2
                evp = IMInternalModes.waveModesAtWavenumber(N2=N2,zDomain=[-1000 0],k=i*1e-4);
                bases{i} = solver.solveEVP(evp,nModes=3);
            end
            bases{2} = bases{2}.addNormalization("double",@(b,j) 2*b.innerProductNormFactor(j));
            bases{2}.normalization = "double";
            collection = IMBasisCollection(bases,basisIndex=[2 1 2]);
        end
    end
end
