classdef TestPortableCompatibilityMatrix < matlab.unittest.TestCase
    properties
        root
        matrix
    end
    methods (TestClassSetup)
        function prepare(testCase)
            testCase.root=string(fileparts(fileparts(mfilename('fullpath'))));
            addpath(fullfile(testCase.root,'tools'));
            testCase.matrix=portableCompatibilityDefinition(testCase.root);
        end
    end
    methods (Test, TestTags={'full'})
        function authoritativeInventoryHasNoUnqualifiedSamplingGaps(testCase)
            result=validatePortableCompatibilityMatrix(testCase.matrix,repositoryRoot=testCase.root,expected=testCase.matrix);
            testCase.verifyGreaterThan(result.rowCount,666);
            testCase.verifyEqual(result.unqualifiedCount,0);
            testCase.verifyTrue(testCase.matrix.readiness.standardParityReady);
            testCase.verifyEqual(testCase.matrix.readiness.decision,"STANDARD-PORTABLE-PARITY");
            testCase.verifyEqual(validatePortableCompatibilityMatrix(testCase.matrix,repositoryRoot=testCase.root,expected=testCase.matrix,requireComplete=true).unqualifiedCount,0);
        end
        function generatedAssemblyMatchesCommittedFiles(testCase)
            output=string(tempname); mkdir(output); cleanup=onCleanup(@()rmdir(output,'s'));
            result=generatePortableCompatibilityMatrix(repositoryRoot=testCase.root,outputRoot=output);
            testCase.verifyTrue(isequal(fileread(result.catalogPath),fileread(fullfile(testCase.root,'PortableRuntime/contracts/portable-compatibility-matrix-v1.json'))),'Regenerate the compatibility JSON from the final source files.');
            testCase.verifyTrue(isequal(fileread(result.documentationPath),fileread(fullfile(testCase.root,'PortableRuntime/COMPATIBILITY.md'))),'Regenerate COMPATIBILITY.md from the final assembly.');
        end
        function removedAndDuplicateRowsAreRejected(testCase)
            candidate=testCase.matrix; candidate.rows(1)=[]; testCase.reject(candidate);
            candidate=testCase.matrix; candidate.rows(end+1)=candidate.rows(1); testCase.reject(candidate);
        end
        function staleAuthorityAndInventoriesAreRejected(testCase)
            candidate=testCase.matrix; candidate.sources(1).sha256=repmat("0",1,64); testCase.reject(candidate);
            candidate=testCase.matrix; candidate.versions(1).version=999; testCase.reject(candidate);
            candidate=testCase.matrix; candidate.inventory.observers(end+1)="WVNewObserver"; testCase.reject(candidate);
            candidate=testCase.matrix; candidate.inventory.variables(end+1)="new_registered_field"; testCase.reject(candidate);
            candidate=testCase.matrix; candidate.inventory.configurations(end+1)="new-transform-aa0"; testCase.reject(candidate);
        end
        function unresolvedAndContradictoryFixturesAreRejected(testCase)
            candidate=testCase.matrix; candidate.rows(1).fixtures="missing-fixture"; testCase.reject(candidate);
            candidate=testCase.matrix; candidate.witnesses(1).symbol="missingMethod"; testCase.reject(candidate);
            candidate=testCase.matrix;
            index=find(arrayfun(@(w)~isempty(w.coverage.rowIds),candidate.witnesses),1);
            candidate.witnesses(index).coverage.rowIds(1)=[]; testCase.reject(candidate);
        end
        function samplingGapsAreRemovedFromTheStandardSlice(testCase)
            rows=testCase.matrix.rows;
            issueRows=rows([rows.issue]==454);
            testCase.verifyEmpty(issueRows);
            testCase.verifyTrue(testCase.matrix.readiness.standardParityReady);
        end
        function supportedSamplingCannotBecomeAGapOrExclusion(testCase)
            candidate=testCase.matrix;
            index=find(string({candidate.rows.slice})=="fields",1);
            testCase.assertNotEmpty(index);
            candidate.rows(index).status="intentional-incompatibility"; testCase.reject(candidate);
            candidate=testCase.matrix; candidate.rows(index).status="unqualified"; candidate.rows(index).issue=454; candidate.rows(index).reason="matlab-supported-sampling-not-implemented"; testCase.reject(candidate);
            candidate=testCase.matrix; candidate.readiness.standardParityReady=false; candidate.readiness.decision="pending-307"; testCase.reject(candidate);
            candidate=testCase.matrix; candidate.readiness.decision="pending-307"; testCase.reject(candidate);
        end
    end
    methods (Access=private)
        function reject(testCase,candidate)
            testCase.verifyError(@()validatePortableCompatibilityMatrix(candidate,repositoryRoot=testCase.root,expected=testCase.matrix),'WaveVortexModel:InvalidCompatibilityMatrix');
        end
    end
end
