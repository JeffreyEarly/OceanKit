classdef IMNativeBasisProtocolTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function configurePath(testCase)
            original = path;
            root = fileparts(fileparts(mfilename("fullpath")));
            addpath(root,fullfile(root,"UnitTestsV2","TestSupport"));
            testCase.addTeardown(@() path(original));
        end
    end

    methods (Test)
        function numericalAndAnalyticalWavesKeepTheManuscriptConvention(testCase)
            D = 1000; k = 2e-4; N0 = .01;
            evp = IMInternalModes.waveModesAtWavenumber(N2=@(z) N0^2*ones(size(z)),zDomain=[-D 0],k=k,surfaceBoundary=IMBoundaryCondition.dirichlet());
            numeric = IMSolverSpectral(nEVP=48).solveEVP(evp,nModes=4);
            exact = IMConstantStratificationSolution(N0=N0,zDomain=[-D 0]).internalModes(evp,nModes=4);
            z = linspace(-D,0,401).';
            for basis = {numeric,exact}
                b = basis{1};
                testCase.verifyEqual(b.F(z),D*b.Gz(z),AbsTol=1e-10);
                testCase.verifyEqual(b.Fz(z),D*b.Gzz(z),AbsTol=1e-10);
                testCase.verifyEqual(b.evaluate(z,variable="G",derivativeOrder=2),b.Gzz(z),AbsTol=1e-12);
                testCase.verifyEqual(trapz(z,(b.F(z).^2+(k*D)^2*b.G(z).^2))/D,ones(1,4),AbsTol=1e-8);
                testCase.verifyEqual(b.F(z,normalization="maxAbsF"),b.withNormalization("maxAbsF").F(z),AbsTol=0);
            end
        end

        function selectionKeepsCustomScalesAndIndependentCopies(testCase)
            evp = IMInternalModes.waveModesAtWavenumber(N2=@(z) 1e-4*ones(size(z)),zDomain=[-1000 0],k=2e-4,surfaceBoundary=IMBoundaryCondition.dirichlet());
            source = IMSolverSpectral(nEVP=32).solveEVP(evp,nModes=4);
            exact = IMConstantStratificationSolution(N0=.01,zDomain=[-1000 0]).internalModes(evp,nModes=4);
            z = linspace(-1000,0,81).';
            for basis = {source,exact}
                b = basis{1}.addNormalization("custom",@(~,j) j+1).withNormalization("custom");
                F = b.F(z);
                selected = b.selectColumns([4 2]);
                testCase.verifyEqual(class(selected),class(b));
                testCase.verifyEqual(selected.F(z),F(:,[4 2]),AbsTol=1e-12);
                testCase.verifyEqual(selected.columnLabels,b.columnLabels([4 2]));
                changed = selected.addNormalization("custom",@(~,~) 2);
                testCase.verifyNotEqual(changed.F(z),selected.F(z));
                testCase.verifyEqual(b.F(z),F,AbsTol=0);
                testCase.verifySize(b.selectColumns([]).F(z),[81 0]);
            end
        end

        function builtInStyleCacheIsLazyAndDetachedOnRuleChanges(testCase)
            evp = IMInternalModes.waveModesAtWavenumber(N2=@(z) 1e-4*ones(size(z)),zDomain=[-1000 0],k=2e-4,surfaceBoundary=IMBoundaryCondition.dirichlet());
            source = IMSolverSpectral(nEVP=24).solveEVP(evp,nModes=3);
            calls = 0;
            b = IMInternalModesBasisTestAccess(source).withVectorRule("counted",@compute);
            b = b.withNormalization("counted");
            testCase.verifyEqual(calls,0);
            testCase.verifyEqual(b.normalizationFactors(),[1 2 3]);
            testCase.verifyEqual(b.normalizationFactors(),[1 2 3]);
            testCase.verifyEqual(calls,1);
            selected = b.selectColumns([3 1]);
            testCase.verifyEqual(selected.normalizationFactors(),[3 1]);
            testCase.verifyEqual(calls,1);
            changed = b.addNormalization("counted",@(~,~) 7);
            testCase.verifyEqual(changed.normalizationFactors(),[7 7 7]);
            testCase.verifyEqual(b.normalizationFactors(),[1 2 3]);
            testCase.verifyEqual(calls,1);
            function values = compute(basis)
                calls = calls+1;
                values = 1:basis.nColumns;
            end
        end

        function customDiagnosticsRequireTheirOwnDerivativeContract(testCase)
            z = linspace(-1,0,17).';
            solver = IMSolverSpectral(nEVP=32);
            for formulation = ["F","G"]
                evp = IMInternalModes(rotationRate=0,name="customDiagnostic", ...
                    formulation=formulation,N2=@(z) ones(size(z)),zDomain=[-1 0],p=1,q=0,r=1, ...
                    FfromGz=@(~,dGdz,~,~) dGdz.^2+1,GfromFz=@(~,dFdz,~,~) dFdz.^2+1);
                basis = solver.solveEVP(evp,nModes=2);
                testCase.verifyEqual(basis.evaluate(z,variable=formulation,derivativeOrder=1),basis.uz(z),AbsTol=0);
                if formulation == "G", diagnostic = "F"; else, diagnostic = "G"; end
                description = basis.capabilities;
                iVariable = find(description.variables == diagnostic,1);
                testCase.verifyFalse(ismember(1,description.derivativeOrders{iVariable}));
                testCase.verifyError(@() basis.evaluate(z,variable=diagnostic,derivativeOrder=1),"IMInternalModesBasis:UnsupportedDerivative");
            end

            evp = IMInternalModes(rotationRate=0,name="customExactDiagnostic", ...
                formulation="G",N2=@(z) ones(size(z)),zDomain=[-1 0],p=1,q=0,r=1);
            exact = IMAnalyticalInternalModesBasis(evp=evp,h=1/pi^2,modeNumber=1,N2=evp.N2, ...
                rawVariableFunction=@(~,variable,z) (variable=="G").*sin(pi*(z+1))+(variable=="F").*(pi^2*cos(pi*(z+1)).^2+1), ...
                rawUzFunction=@(~,z) pi*cos(pi*(z+1)));
            testCase.verifyEqual(exact.Gz(z),exact.uz(z),AbsTol=0);
            description = exact.capabilities;
            testCase.verifyFalse(ismember(1,description.derivativeOrders{description.variables == "F"}));
            testCase.verifyFalse(ismember(2,description.derivativeOrders{description.variables == "G"}));
            testCase.verifyError(@() exact.Fz(z),"IMAnalyticalInternalModesBasis:UnsupportedDerivative");
            testCase.verifyError(@() exact.Gzz(z),"IMAnalyticalInternalModesBasis:UnsupportedDerivative");
        end

        function storedCompanionDerivativesStayAlignedWithEvaluation(testCase)
            N2 = @(z) 1e-4*ones(size(z));
            evp = IMInternalModes.geostrophicAPVModes(rotationRate=0,N2=N2,zDomain=[-1000 0],g0=0.02,gd=Inf,surfaceBoundary="freeSurface");
            z = linspace(-1000,0,33).';
            for solver = {IMSolverSpectral(nEVP=48),IMSolverFiniteDifference(z=linspace(-1000,0,65).')}
                basis = solver{1}.solveEVP(evp,nModes=3);
                testCase.verifyEqual(basis.evaluate(z,variable="G",derivativeOrder=1),basis.Gz(z),AbsTol=0);
            end

            mda = IMSolverSpectral(nEVP=48).solveEVP(IMInternalModes.meanDensityAnomalyModes( ...
                rotationRate=0,N2=N2,zDomain=[-1000 0],g0=Inf,gd=Inf),nModes=3);
            expectedFz = -N2(z)/mda.evp.g.*mda.G(z);
            testCase.verifyEqual(mda.Fz(z),expectedFz,AbsTol=1e-10);
            testCase.verifyEqual(mda.evaluate(z,variable="F",derivativeOrder=1),mda.Fz(z),AbsTol=0);
        end
    end
end
