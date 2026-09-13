classdef TestPortableHistoricalQualification < matlab.unittest.TestCase
    properties (SetAccess=private)
        root (1,1) string
        inventory (1,1) struct
    end
    methods (TestClassSetup)
        function loadAuthority(testCase)
            testCase.root = string(fileparts(fileparts(mfilename("fullpath"))));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(testCase.root,"tools")));
            testCase.inventory = jsondecode(fileread(fullfile(testCase.root,"PortableRuntime","qualification","historical-test-inventory-v1.json")));
        end
    end
    methods (Test,TestTags="full")
        function historicalWorkloadSurvivesCurrentClassGrowth(testCase)
            for entry = reshape(testCase.inventory.entries,1,[])
                [report,validator,identifier] = testCase.family(entry);
                currentNames = strings(1,0);
                for source = reshape(entry.sources,1,[])
                    suite = testsuite(fullfile(testCase.root,source.path));
                    currentNames = [currentNames,string({suite.Name})]; %#ok<AGROW>
                end
                testCase.assertTrue(any(~ismember(currentNames,string(entry.testNames))),"Fixture requires actual test-class growth since original execution.");
                validation = validator(report,evidenceScope="historical",repositoryRoot=testCase.root);
                testCase.verifyEqual(validation.scope,"historical-workload");
                testCase.verifyFalse(validation.currentReadiness);
                testCase.verifyEqual(validation.sourceCommit,string(entry.sourceCommit));
                testCase.verifyError(@()validator(report,repositoryRoot=testCase.root),identifier);
            end
        end
        function missingDuplicateOrFailedOriginalTestsAreRejected(testCase)
            for entry = reshape(testCase.inventory.entries,1,[])
                [report,validator,identifier] = testCase.family(entry);
                bad = report; bad.tests(end) = [];
                testCase.verifyError(@()validator(bad,evidenceScope="historical"),identifier);
                bad = report; bad.tests(end+1) = bad.tests(1);
                testCase.verifyError(@()validator(bad,evidenceScope="historical"),identifier);
                for flag = ["passed","failed","incomplete"]
                    bad = report; bad.tests(1).(flag) = flag~="passed";
                    testCase.verifyError(@()validator(bad,evidenceScope="historical"),identifier);
                end
                bad = report; bad.elapsedSeconds = bad.elapsedSeconds+1;
                testCase.verifyError(@()validator(bad,evidenceScope="historical"),identifier);
            end
        end
        function historicalResultsCannotBePromotedByAddingCurrentTests(testCase)
            for entry = reshape(testCase.inventory.entries,1,[])
                [report,validator,identifier] = testCase.family(entry);
                synthetic = report;
                synthetic.sourceCommit = repmat('a',1,40);
                [currentNames,current] = portableQualificationTestInventory(synthetic,string(entry.family),testCase.root,"current");
                testCase.verifyEqual(current.scope,"current-inventory-only");
                testCase.verifyFalse(current.currentReadiness);
                testCase.verifyGreaterThan(numel(currentNames),numel(report.tests));
                % Even an otherwise current report must include the grown
                % inventory. This structural fixture is not execution evidence.
                testCase.verifyError(@()validator(synthetic),identifier);
                promoted = report;
                for missing = reshape(setdiff(currentNames,string({report.tests.name})),1,[])
                    promoted.tests(end+1) = report.tests(1);
                    promoted.tests(end).name = char(missing);
                end
                testCase.verifyError(@()validator(promoted),identifier);
                testCase.verifyError(@()validator(promoted,evidenceScope="historical"),identifier);
                % Relabeling a synthetic complete inventory with an unknown
                % commit cannot establish actual source/execution freshness.
                promoted.sourceCommit = repmat('a',1,40);
                inventoryOnly = validator(promoted);
                testCase.verifyEqual(inventoryOnly.scope,"current-inventory-only");
                testCase.verifyFalse(inventoryOnly.currentReadiness);
                promoted = report; promoted.sourceCommit = repmat('a',1,40);
                testCase.verifyError(@()validator(promoted,evidenceScope="historical"),identifier);
            end
        end
        function committedAuthorityRejectsTamperedManifestAndArtifact(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            fixtureRoot = string(fixture.Folder);
            target = fullfile(fixtureRoot,"PortableRuntime","qualification");
            mkdir(target);
            inventoryPath = "PortableRuntime/qualification/historical-test-inventory-v1.json";
            copyfile(fullfile(testCase.root,inventoryPath),fullfile(fixtureRoot,inventoryPath));
            entry = testCase.inventory.entries(1);
            [report,~,identifier] = testCase.family(entry);
            copyfile(fullfile(testCase.root,entry.reportPath),fullfile(fixtureRoot,entry.reportPath));
            [names,result] = portableQualificationTestInventory(report,string(entry.family),fixtureRoot,"historical");
            testCase.verifyEqual(sort(names),sort(reshape(string(entry.testNames),1,[])));
            testCase.verifyFalse(result.currentReadiness);
            changedInventory = testCase.inventory;
            changedInventory.entries(1).testNames(end) = [];
            writeJSON(fullfile(fixtureRoot,inventoryPath),changedInventory);
            testCase.verifyError(@()portableQualificationTestInventory(report,string(entry.family),fixtureRoot,"historical"),identifier);
            copyfile(fullfile(testCase.root,inventoryPath),fullfile(fixtureRoot,inventoryPath));
            changed = report; changed.tests(1).passed = false;
            writeJSON(fullfile(fixtureRoot,entry.reportPath),changed);
            testCase.verifyError(@()portableQualificationTestInventory(report,string(entry.family),fixtureRoot,"historical"),identifier);
            delete(fullfile(fixtureRoot,entry.reportPath));
            testCase.verifyError(@()portableQualificationTestInventory(report,string(entry.family),fixtureRoot,"historical"),identifier);
        end
    end
    methods (Access=private)
        function [report,validator,identifier] = family(testCase,entry)
            report = jsondecode(fileread(fullfile(testCase.root,entry.reportPath)));
            switch string(entry.family)
                case "stratified-qg"
                    validator = @validatePortableStratifiedQGQualification;
                    identifier = "WaveVortexModel:InvalidSQGQualification";
                case "hydrostatic"
                    validator = @validatePortableHydrostaticQualification;
                    identifier = "WaveVortexModel:InvalidHydrostaticQualification";
                case "boussinesq"
                    validator = @validatePortableBoussinesqQualification;
                    identifier = "WaveVortexModel:InvalidBoussinesqQualification";
            end
        end
    end
end

function writeJSON(path,value)
stream = fopen(path,'w');
assert(stream>=0,"Cannot open historical qualification test artifact.");
cleanup = onCleanup(@()fclose(stream));
fprintf(stream,"%s\n",jsonencode(value,PrettyPrint=true));
end
