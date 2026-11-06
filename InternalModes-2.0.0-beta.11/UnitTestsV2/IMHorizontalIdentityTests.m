classdef IMHorizontalIdentityTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function configurePath(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,"UnitTestsV2","TestSupport")));
        end
    end
    methods (Test)
        function fixedFrequencyRootsSatisfyRigidLidDispersion(testCase)
            N0 = .01; D = 1000;
            for latitude = [-33 33]
                bg = IMBackground.constantStratification(N0=N0,zDomain=[-D 0],latitude=latitude,surfaceCondition="rigidLid");
                for omega = [.5*abs(bg.f0),.005]
                    for solver = {IMSolverAnalytical(),IMSolverSpectral(nEVP=48)}
                        b = createIMTestBasis("frequency",bg,omega=omega,solver=solver{1},nModes=3);
                        expected = (omega^2-bg.f0^2)*(pi*(1:3)/D).^2/(N0^2-omega^2);
                        testCase.verifyEqual(b.horizontalKind,"perColumn");
                        testCase.verifyEmpty(b.k);
                        testCase.verifyEqual(b.kSquared,expected,RelTol=3e-11);
                    end
                end
            end
        end
        function selectingRootsPreservesDispersionAndScientificIdentity(testCase)
            bg = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            b = createIMTestBasis("frequency",bg,nModes=4,solver=IMSolverAnalytical()); s = b.selectColumns([4 2]);
            testCase.verifyEqual(s.kSquared,b.kSquared([4 2]),AbsTol=0);
            testCase.verifyEqual(s.scientificSpecification.resolved.horizontal.kSquared,s.kSquared);
            testCase.verifyEqual(s.constructionRecord.requested,b.constructionRecord.requested);
            testCase.verifyEqual(s.selectColumns([]).horizontalKind,"perColumn");
        end
        function independentFamiliesAcceptExternallyAssignedPages(testCase)
            bg = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            for family = ["APV","MDA","hydrostaticF","hydrostaticG"]
                b = createIMTestBasis(family,bg);
                testCase.verifyEqual(b.horizontalKind,"independent"); testCase.verifyEmpty(b.k);
                c = IMBasisCollection({b},kappa=[1e-3 3e-3 1e-3],basisIndex=[1 1 1]);
                testCase.verifyEqual(c.kappa,[1e-3 3e-3 1e-3]);
                testCase.verifyEqual(c.metadata.horizontalIdentity.kind,"independent");
            end
        end
        function solvedWavenumbersCannotBeReassigned(testCase)
            bg = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            for family = ["wave","inertial","enstrophy"]
                b = createIMTestBasis(family,bg);
                testCase.verifyEqual(b.horizontalKind,"fixedWavenumber");
                testCase.verifyEqual(IMBasisCollection({b}).kappa,b.k);
                testCase.verifyError(@() IMBasisCollection({b},kappa=b.k+1e-3),"IMBasisCollection:KappaMismatch");
            end
            inertial = createIMTestBasis("inertial",bg); testCase.verifyEqual(inertial.k,0);
        end
        function perColumnRootsRemainDistinctFromCollectionPages(testCase)
            bg = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            b = createIMTestBasis("frequency",bg,solver=IMSolverAnalytical());
            c = IMBasisCollection({b},basisIndex=[1 1]); testCase.verifyEmpty(c.kappa);
            testCase.verifyEqual(c.metadata.horizontalIdentity.kSquared,b.kSquared);
            testCase.verifyError(@() IMBasisCollection({b},kappa=1e-3),"IMBasisCollection:UnsupportedHorizontalPage");
            empty = b.selectColumns([]);
            testCase.verifyError(@() IMBasisCollection({empty},kappa=1e-3),"IMBasisCollection:UnsupportedHorizontalPage");
        end
        function developerScalarBasesDeclareUnspecifiedHorizontalMeaning(testCase)
            problem = IMEigenvalueProblem(zDomain=[-1 0],p=1,q=0,r=1,surfaceBoundary=IMBoundaryCondition.dirichlet(),bottomBoundary=IMBoundaryCondition.dirichlet());
            b = IMSolverSpectral(nEVP=32).solveEVP(problem,nModes=2);
            testCase.verifyEqual(b.horizontalKind,"unspecified"); testCase.verifyEmpty(b.k);
            testCase.verifyEqual(b.capabilities.variables,"u");
            testCase.verifyError(@() IMBasisCollection({b},kappa=1e-3),"IMBasisCollection:UnsupportedHorizontalPage");
        end
    end
end
