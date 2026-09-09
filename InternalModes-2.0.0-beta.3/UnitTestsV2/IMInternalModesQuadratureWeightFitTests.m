classdef IMInternalModesQuadratureWeightFitTests < matlab.unittest.TestCase

    properties (Access = private)
        originalPath
        basisSet
        z
    end

    methods (TestClassSetup)
        function configureFamily(testCase)
            repoRoot = fileparts(fileparts(mfilename("fullpath")));
            testCase.originalPath = path;
            addpath(repoRoot);
            D = 1200;
            N2 = @(z) 1e-4*exp(z/2000);
            solver = IMSolverSpectral(nEVP=128,coordinateKind="wkb");
            testCase.basisSet = solver.solveEVP(IMInternalModes.hydrostaticFModes(N2=N2,zDomain=[-D 0]),nModes=8);
            sigma = linspace(0,1,19).';
            testCase.z = -D+D*(1-(1-sigma).^2);
        end
    end

    methods (TestClassTeardown)
        function restorePath(testCase)
            path(testCase.originalPath);
        end
    end

    methods (Test)
        function stackedSystemMatchesIndependentVariableRows(testCase)
            [weights,fit] = testCase.basisSet.quadratureWeightsForPoints(z=testCase.z,nModes=5,variables=["F","G"]);
            testCase.verifyClass(fit,"IMInternalModesQuadratureWeightFit")
            testCase.verifyEqual(fit.availableVariables,["F","G"])
            testCase.verifyEqual(weights,fit.transform.weights,AbsTol=0)
            testCase.verifyEqual(sum(weights),1200,AbsTol=1e-8)
            testCase.verifyGreaterThanOrEqual(min(weights),-1e-10)
            testCase.verifyEqual(fit.residualNorm,norm(fit.residual,2),AbsTol=1e-13)
            testCase.verifyEqual(fit.geometricResidualNorm,norm(fit.geometricResidual,2),AbsTol=1e-13)
            testCase.verifyEqual(fit.variableResidualNorm(variable="F"),norm(fit.residual(fit.objectiveRowVariables == "F")),AbsTol=1e-13)
            testCase.verifyEqual(fit.variableResidualNorm(variable="G"),norm(fit.residual(fit.objectiveRowVariables == "G")),AbsTol=1e-13)
            testCase.verifyEqual(sum(fit.objectiveRowVariables == "F"),15)
            testCase.verifyEqual(sum(fit.objectiveRowVariables == "G"),10)
            testCase.verifyLessThanOrEqual(fit.residualNorm,fit.geometricResidualNorm+1e-10)
        end

        function objectiveRowsReconstructBothNormalizedGramSystems(testCase)
            [~,fit] = testCase.basisSet.quadratureWeightsForPoints(z=testCase.z,nModes=4,variables=["F","G"]);
            for variable = ["F","G"]
                rows = fit.objectiveRowVariables == variable;
                transform = fit.transform;
                active = transform.activeModeMask(variable=variable);
                target = transform.targetGramMatrix(variable=variable);
                gram = transform.gramMatrix(variable=variable);
                scale = 1./sqrt(abs(diag(target(active,active))));
                expected = norm(scale.*(gram(active,active)-target(active,active)).*scale.',"fro");
                testCase.verifyEqual(norm(fit.residual(rows)),expected,RelTol=1e-11,AbsTol=1e-12)
            end
        end

        function customObjectiveReceivesStackAndPerVariableContexts(testCase)
            objective = @(context) testCase.validateAndRegularize(context);
            [~,fit] = testCase.basisSet.quadratureWeightsForPoints(z=testCase.z,nModes=4,variables=["F","G"],objective=objective);
            testCase.verifyEqual(fit.objectiveName,"stackedWithTinyGeometricPenalty")
            testCase.verifyTrue(any(fit.objectiveRowVariables == "custom"))
            testCase.verifyEqual(size(fit.objectiveMatrix,2),length(testCase.z))
        end

        function explicitWeightsBypassSharedFitAndPreserveProvenance(testCase)
            geometric = geometricWeights(testCase.z,testCase.basisSet.zDomain);
            [transform,assessment] = testCase.basisSet.discreteTransform(z=testCase.z,weights=geometric,nModes=4,variables=["F","G"],gramTolerance=10);
            testCase.verifyEmpty(assessment.weightFit)
            testCase.verifyEqual(transform.weights,geometric,AbsTol=0)
            testCase.verifyEqual(assessment.requestedPointCount,length(testCase.z))
            testCase.verifyEqual(assessment.actualPointCount,length(testCase.z))
            testCase.verifyEqual(assessment.availableVariables,["F","G"])
        end

        function zeroActiveChannelUsesGeometricRuleAndKeepsFamilyAlignment(testCase)
            [weights,fit] = testCase.basisSet.quadratureWeightsForPoints(z=testCase.z,nModes=1,variables="G");
            expected = geometricWeights(testCase.z,testCase.basisSet.zDomain);
            testCase.verifyEqual(weights,expected,AbsTol=0)
            testCase.verifyEmpty(fit.objectiveMatrix)
            testCase.verifyEqual(fit.residualNorm,0,AbsTol=0)
            testCase.verifyEqual(fit.transform.activeModeMask(variable="G"),false)
            testCase.verifyEqual(fit.transform.forwardMatrix(variable="G"),zeros(1,length(testCase.z)),AbsTol=0)
            values = fit.transform.transformBack(3,variable="G");
            testCase.verifyEqual(values,zeros(size(testCase.z)),AbsTol=1e-12)
        end
    end

    methods (Access = private)
        function specification = validateAndRegularize(testCase,context)
            testCase.verifyEqual(context.availableVariables,["F","G"])
            testCase.verifyTrue(isfield(context.variableContexts,"F"))
            testCase.verifyTrue(isfield(context.variableContexts,"G"))
            testCase.verifyEqual(length(context.normalizedGramVariables),size(context.normalizedGramA,1))
            lambda = 1e-12;
            penalty = sqrt(lambda)*diag(1./context.geometricWeights);
            specification = struct("A",[context.normalizedGramA;penalty], ...
                "b",[context.normalizedGramB;sqrt(lambda)*ones(length(context.z),1)], ...
                "name","stackedWithTinyGeometricPenalty", ...
                "rowVariables",[context.normalizedGramVariables;repmat("custom",length(context.z),1)], ...
                "modePairs",[context.normalizedGramModePairs;ones(length(context.z),2)]);
        end
    end
end

function weights = geometricWeights(z,zDomain)
edges = [zDomain(1);0.5*(z(1:end-1)+z(2:end));zDomain(2)];
weights = diff(edges);
end
