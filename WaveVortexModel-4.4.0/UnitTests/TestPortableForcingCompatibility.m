classdef TestPortableForcingCompatibility < matlab.unittest.TestCase
    properties (SetAccess=private)
        root (1,1) string
        matrix (1,1) struct
    end
    methods (TestClassSetup)
        function loadCatalog(testCase)
            testCase.root = string(fileparts(fileparts(mfilename("fullpath"))));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(testCase.root,"tools")));
            testCase.matrix = jsondecode(fileread(fullfile(testCase.root,"PortableRuntime","contracts","portable-forcing-compatibility-v1.json")));
        end
    end
    methods (Test,TestTags="full")
        function matlabRegenerationMatchesCommittedArtifacts(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            result = generatePortableForcingCompatibility(outputRoot=string(fixture.Folder));
            testCase.verifyEqual(fileread(result.catalogPath),fileread(fullfile(testCase.root,"PortableRuntime","contracts","portable-forcing-compatibility-v1.json")));
            testCase.verifyEqual(fileread(result.headerPath),fileread(fullfile(testCase.root,"PortableRuntime","tests","generated","WVForcingCompatibilityRows.hpp")));
        end
        function incompleteEvidenceCannotClaimReadiness(testCase)
            validatePortableForcingCompatibility(testCase.matrix,requireComplete=true);
            pending = testCase.matrix;
            pending.rows(1).acceptance = 'pending';
            pending.completion = 'incomplete';
            testCase.verifyError(@()validatePortableForcingCompatibility(pending,requireComplete=true),"WaveVortexModel:CompatibilityIncomplete");
            pending.completion = 'complete';
            testCase.verifyInvalid(pending);
            bad = testCase.matrix;
            bad.rows(1).evidence = {'matlab-contracts'};
            testCase.verifyInvalid(bad);
            bad = testCase.matrix;
            bad.rows(1).acceptance = 'incompatible';
            testCase.verifyInvalid(bad);
        end
        function invalidIdentitiesAndEvidenceAreRejected(testCase)
            bad = testCase.matrix;
            bad.rows(end) = [];
            testCase.verifyInvalid(bad);
            bad = testCase.matrix;
            bad.rows(end) = bad.rows(1);
            testCase.verifyInvalid(bad);
            bad = testCase.matrix;
            bad.rows(1).contractVersion = 2;
            testCase.verifyInvalid(bad);
            bad = testCase.matrix;
            bad.rows(1).evidence = {'does-not-exist'};
            testCase.verifyInvalid(bad);
            bad = testCase.matrix;
            bad.evidence(1).symbol = 'doesNotExist';
            testCase.verifyInvalid(bad);
            bad = testCase.matrix;
            bad.evidence(1).path = '../outside.m';
            testCase.verifyInvalid(bad);
            bad = testCase.matrix;
            bad.configurations(1).shouldAntialias = ~bad.configurations(1).shouldAntialias;
            testCase.verifyInvalid(bad);
        end
        function implementedPairsAndScientificRejectionsStayDistinct(testCase)
            rows = testCase.matrix.rows;
            for family = ["constant-hydrostatic","constant-nonhydrostatic","barotropic","stratified-qg","hydrostatic","boussinesq"]
                antialias = testCase.row(family+"-aa0/WVAntialiasing");
                testCase.verifyEqual(string(antialias.matlab.applicability),"applicable");
                testCase.verifyEqual(antialias.matlab.priority,127);
                testCase.verifyTrue(antialias.implementation.factoryAvailable);
                testCase.verifyEqual(antialias.implementation.issue,0);
                doubleAntialias = testCase.row(family+"-aa1/WVAntialiasing");
                testCase.verifyEqual(string(doubleAntialias.matlab.rejection),"double-antialias");
            end
            barotropic = testCase.row("barotropic-aa0/WVVerticalDiffusivity");
            testCase.verifyEqual(string(barotropic.matlab.rejection),"forcing-stage-unavailable");
            wave = testCase.row("constant-hydrostatic-aa0/WVVerticalDiffusivity");
            testCase.verifyEqual(string(wave.matlab.applicability),"applicable");
            testCase.verifyEqual(wave.implementation.issue,0);
            testCase.verifyEqual(numel(rows),144);
            testCase.verifyEqual(string(testCase.matrix.inventory.excluded.identity),"WVThermalDamping");
            testCase.verifyFalse(any(string({rows.forcing}) == "WVThermalDamping"));
        end
    end
    methods (Access=private)
        function row = row(testCase,id)
            row = testCase.matrix.rows(string({testCase.matrix.rows.id}) == id);
            testCase.assertTrue(isscalar(row));
        end
        function verifyInvalid(testCase,matrix)
            testCase.verifyError(@()validatePortableForcingCompatibility(matrix),"WaveVortexModel:InvalidForcingCompatibility");
        end
    end
end
