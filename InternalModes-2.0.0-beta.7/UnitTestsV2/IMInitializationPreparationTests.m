classdef IMInitializationPreparationTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addSupport(testCase)
            folder=fileparts(mfilename('fullpath'));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(folder,'TestSupport')));
        end
    end
    methods (Test)
        function vectorizedReportsMatchFrozenScalarReference(testCase)
            z=linspace(-1000,0,101).'; w=1+(1:101).'/101;
            x=z/1000; modes=1:12;
            identity=struct(family="waves",columnLabels=string(modes),normalization="test",zDomain=[-1000 0],kappa=1e-4);
            a=struct(identity=identity,values=struct(F=cos(x*modes),G=sin(x*modes)),derivatives=struct(F=-sin(x*modes).*modes/1000,G=cos(x*modes).*modes/1000),equivalentDepths=1./modes,provenance=struct(nEVP=65));
            for scenario=1:9
                c=a; b=a;
                switch scenario
                    case 1
                        b.values.G=b.values.G+1e-8*cos(x*modes);
                    case 2
                        permutation=12:-1:1;
                        b.identity.columnLabels=b.identity.columnLabels(permutation);
                        b.equivalentDepths=b.equivalentDepths(permutation);
                        for group=["values","derivatives"]
                            for field=["F","G"], b.(group).(field)=-b.(group).(field)(:,permutation); end
                        end
                    case 3
                        c.identity.columnLabels(2)=c.identity.columnLabels(1);
                        b.identity.columnLabels(4)=b.identity.columnLabels(3);
                        b.identity.columnLabels(5)="missing";
                    case 4
                        c.values.F(:)=0; c.values.G(:)=0; b.values=c.values;
                        b.derivatives.F=-c.derivatives.F; b.derivatives.G=-c.derivatives.G;
                    case 5
                        b.derivatives=rmfield(b.derivatives,'F');
                        b=rmfield(b,'equivalentDepths');
                    case 6
                        c.equivalentDepths=[0 1 Inf -Inf NaN Inf 1 0 2 3 4 5];
                        b.equivalentDepths=[0 0 Inf -Inf 1 -Inf Inf NaN 2 4 4 5];
                    case 7
                        c.values.F(:)=1; c.values.G(:)=1;
                        b.values.F=c.values.F; b.values.G=-c.values.G;
                    case 8
                        c.identity.columnLabels(:)="duplicate";
                    case 9
                        c.identity.zDomain=[-realmax realmax]; b.identity=c.identity;
                        c.derivatives.F(:)=0; b.derivatives.F(:)=0;
                end
                expected=assessModeConvergenceScalarReference(c,b,z,w);
                actual=assessModeConvergence(c,b,z,w);
                testCase.verifyEqual(rmfield(actual,'costs'),rmfield(expected,'costs'),AbsTol=1e-13);
            end
        end

        function preparedDerivativesMatchIndependentEvaluation(testCase)
            for kind=["z","wkb","density"]
                firstEVP=IMInternalModes.hydrostaticGModes(N2=@(z)1e-4*exp(2*z/700),zDomain=[-1000 0]);
                secondEVP=IMInternalModes.hydrostaticGModes(N2=@(z)4e-4+0*z,zDomain=[-500 0]);
                first=IMSolverSpectral(nEVP=33,coordinateKind=kind).configuredForEVP(firstEVP);
                original=first.differentiateGridValues(eye(33),2);
                second=first.configuredForEVP(secondEVP);
                fresh=IMSolverSpectral(nEVP=33,coordinateKind=kind).configuredForEVP(secondEVP);
                testCase.verifyEqual(second.differentiateGridValues(eye(33),2),fresh.differentiateGridValues(eye(33),2));
                testCase.verifyEqual(first.differentiateGridValues(eye(33),2),original);
                for solver={first,second}
                    s=solver{1}; values=cos((1:33).'*(1:7)/33);
                    for order=0:2
                        expected=values;
                        if order>0
                            expected=s.evaluatePhysicalDerivative(s.T\values,s.zNative,order);
                        end
                        testCase.verifyEqual(s.differentiateGridValues(values,order),expected,AbsTol=1e-11,RelTol=1e-11);
                    end
                end
            end
        end

        function repeatedCallsReusePreparationAndPreserveSubclassFallback(testCase)
            evp=IMInternalModes.waveModesAtWavenumber(N2=@(z)1e-4+0*z,zDomain=[-1000 0],k=1e-4,surfaceBoundary=IMBoundaryCondition(a=0,b=1,c=1,d=0));
            solver=IMSolverSpectral(nEVP=33).configuredForEVP(evp);
            basis=solver.solveEVP(evp,nModes=4);
            custom=IMBulkUnsupportedSpectralSolver(nEVP=33).solveEVP(evp,nModes=4);
            profile clear
            profile on
            cleanup=onCleanup(@() profile('off'));
            for iteration=1:3
                solver.differentiateGridValues(eye(33),1);
                evp.endpointWeights();
            end
            info=profile('info');
            profile off
            functions=info.FunctionTable;
            names=string({functions.FunctionName});
            testCase.verifyFalse(any(endsWith(names,'>decomposition.decomposition') | names=="decomposition.decomposition"));
            testCase.verifyFalse(any(contains(names,'isEigenvalueDependent')));
            profile clear
            profile on
            actual=basis.endpointGramTerms(useNormalized=false);
            info=profile('info');
            profile off
            testCase.verifyFalse(any(contains(string({info.FunctionTable.FunctionName}),'rawUz')));
            profile clear
            profile on
            expected=custom.endpointGramTerms(useNormalized=false);
            info=profile('info');
            profile off
            clear cleanup
            testCase.verifyEqual(actual,expected,AbsTol=1e-12);
            testCase.verifyTrue(any(contains(string({info.FunctionTable.FunctionName}),'rawUz')));
        end

        function endpointRecipesAndTermsMatchDirectCalculation(testCase)
            boundaries={IMBoundaryCondition.dirichlet(),IMBoundaryCondition(a=0,b=1,c=1,d=0),IMBoundaryCondition(a=1,b=2,c=3,d=1),IMBoundaryCondition(a=1,b=1,c=1,d=1)};
            for index=1:numel(boundaries)
                boundary=boundaries{index};
                evp=IMEigenvalueProblem(surfaceBoundary=boundary,bottomBoundary=boundary);
                for location=["surface","bottom"]
                    terms=evp.endpointWeights(location);
                    coefficient=boundary.endpointWeightCoefficient(location);
                    if boundary.isEigenvalueDependent() && isfinite(coefficient)
                        testCase.verifyEqual(terms,struct(location=location,coefficient=coefficient,c=boundary.c,d=boundary.d));
                    else
                        testCase.verifyEmpty(terms);
                    end
                end
            end
            for boundary=boundaries(2:3)
                evp=IMInternalModes.waveModesAtWavenumber(N2=@(z)1e-4+0*z,zDomain=[-1000 0],k=1e-4,surfaceBoundary=boundary{1});
                basis=IMSolverSpectral(nEVP=33).solveEVP(evp,nModes=4);
                for normalized=[false true]
                    terms=basis.endpointGramTerms(useNormalized=normalized);
                    factors=ones(1,4);
                    if normalized, factors=basis.normalizationFactors(basis.normalization); end
                    weights=evp.endpointWeights("surface");
                    context=evp.contextForSolver(basis.solver);
                    p=IMEigenvalueProblem.evaluateCoefficient(evp.p,0,context);
                    expected=weights.c*(basis.rawU(0)./factors)-weights.d*p*(basis.rawUz(0)./factors);
                    testCase.verifyEqual(terms(1).values,expected,AbsTol=1e-12);
                    testCase.verifyEmpty(basis.endpointGramTerms(zBounds=[-900 -100],useNormalized=normalized));
                end
            end
        end
    end
end
