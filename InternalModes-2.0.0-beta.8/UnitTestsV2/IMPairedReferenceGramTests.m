classdef IMPairedReferenceGramTests < matlab.unittest.TestCase
    % Check physical forms against integration of the evaluated F/G profiles.
    methods (TestClassSetup)
        function configurePath(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
        end
    end

    methods (Test)
        function nearCapacityHydrostaticFormsResolveSpectralProducts(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            basis = IMSolverSpectral(nEVP=24).solveEVP(IMInternalModes.hydrostaticGModes(background),nModes=22);
            [z,weights] = legpts(96,background.zDomain);
            verifyHydrostaticForms(testCase,basis,z,weights(:),background);
        end

        function finiteDifferenceFormsIntegrateReconstructedProfiles(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            knots = linspace(-1000,0,24).';
            basis = IMSolverFiniteDifference(z=knots).solveEVP(IMInternalModes.hydrostaticGModes(background),nModes=22);
            [nodes,rule] = legpts(4,[0 1]);
            z = reshape(knots(1:end-1).'+nodes(:).*diff(knots).',[],1);
            weights = reshape(rule(:).*diff(knots).',[],1);
            verifyHydrostaticForms(testCase,basis,z,weights,background);
        end
    end
end

function verifyHydrostaticForms(testCase,basis,z,weights,background)
for variable = ["F","G"]
    values = basis.evaluate(z,variable=variable);
    if variable == "F", weight = 1; else, weight = .01^2/background.g; end
    expected = values.'*((weight*weights).*values);
    scale = diag(1./sqrt(diag(expected)));
    testCase.verifyEqual(scale*basis.gramMatrix(variable=variable)*scale,scale*expected*scale,AbsTol=2e-11);
    testCase.verifyEqual(scale*basis.majorantGramMatrix(variable=variable)*scale,scale*expected*scale,AbsTol=2e-11);
    physical = basis.projectionOnGrid(z,weights,variable=variable);
    testCase.verifyEqual(physical.projectPhysical(values),eye(basis.nColumns),AbsTol=2e-10);
    columns = [22 19 3];
    selected = basis.selectColumns(columns).withNormalization("rmsG");
    selectedValues = selected.evaluate(z,variable=variable);
    selectedPhysical = selected.projectionOnGrid(z,weights,variable=variable);
    testCase.verifyEqual(selectedPhysical.projectPhysical(selectedValues),eye(3),AbsTol=2e-11);
end
end
